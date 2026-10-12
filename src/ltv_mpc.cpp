#include "bpt/control.hpp"

#include "bpt/bicycle.hpp"
#include "bpt/qp_solver.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <vector>

namespace bpt {
namespace {

Input rate_limited(const Input& desired, const Input& previous, const VehicleParams& params, double dt) {
    Input out;
    const double d_lo = std::max(-params.delta_max, previous.delta - params.delta_rate_max * dt);
    const double d_hi = std::min(params.delta_max, previous.delta + params.delta_rate_max * dt);
    const double a_lo = std::max(params.accel_min, previous.accel - params.accel_rate_max * dt);
    const double a_hi = std::min(params.accel_max, previous.accel + params.accel_rate_max * dt);
    out.delta = std::clamp(desired.delta, d_lo, d_hi);
    out.accel = std::clamp(desired.accel, a_lo, a_hi);
    return out;
}

std::array<double, 6> state_error(const VehicleState& state, const VehicleState& reference) {
    std::array<double, 6> error{};
    const std::array<double, 6> x = state.to_array();
    const std::array<double, 6> r = reference.to_array();
    for (int i = 0; i < 6; ++i) {
        error[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)] - r[static_cast<std::size_t>(i)];
    }
    error[2] = wrap_to_pi(state.psi - reference.psi);
    return error;
}

}  // namespace

LtvMpc::LtvMpc(const Path& path, VehicleParams vehicle, MpcParams tuning, double speed_ref)
    : path_(path), vehicle_(vehicle), tuning_(tuning), speed_ref_(speed_ref) {}

void LtvMpc::reset(const VehicleState& state) {
    s_hint_ = -1.0;
    previous_.delta = state.delta;
    previous_.accel = state.accel;
    have_warm_ = false;
    last_solution_.clear();
}

Command LtvMpc::step(const VehicleState& measured) {
    const auto t0 = std::chrono::steady_clock::now();
    const int horizon = std::max(2, tuning_.horizon);
    const int nvar = 2 * horizon;
    const double dt = tuning_.dt;

    const Path::Projection proj = path_.project(measured.x, measured.y, s_hint_, 30.0);
    s_hint_ = proj.s;

    std::vector<VehicleState> reference(static_cast<std::size_t>(horizon + 1));
    std::vector<Input> reference_input(static_cast<std::size_t>(horizon));
    for (int k = 0; k <= horizon; ++k) {
        const double s = proj.s + static_cast<double>(k) * speed_ref_ * dt;
        const PathPoint point = path_.sample(s);
        const double delta_ss =
            std::clamp(std::atan(vehicle_.wheelbase * point.kappa), -vehicle_.delta_max, vehicle_.delta_max);
        VehicleState ref;
        ref.x = point.x;
        ref.y = point.y;
        ref.psi = point.psi;
        ref.v = speed_ref_;
        ref.delta = delta_ss;
        ref.accel = 0.0;
        reference[static_cast<std::size_t>(k)] = ref;
        if (k < horizon) {
            reference_input[static_cast<std::size_t>(k)] = Input{delta_ss, 0.0};
        }
    }

    // e = S u + w, with u the stacked steering/acceleration commands.
    Matrix sensitivity(6, nvar, 0.0);
    std::vector<double> affine(6);
    const std::array<double, 6> initial_error = state_error(measured, reference[0]);
    for (int i = 0; i < 6; ++i) {
        affine[static_cast<std::size_t>(i)] = initial_error[static_cast<std::size_t>(i)];
    }

    Matrix hessian(nvar, nvar, 0.0);
    std::vector<double> gradient(static_cast<std::size_t>(nvar), 0.0);
    const double stage_weight[4] = {tuning_.w_cte, tuning_.w_along, tuning_.w_heading, tuning_.w_speed};

    auto accumulate = [&](const Matrix& map, const std::vector<double>& offset, double psi_ref, double scale) {
        const double sine = std::sin(psi_ref);
        const double cosine = std::cos(psi_ref);
        double output_map[4][6] = {};
        output_map[0][0] = -sine;
        output_map[0][1] = cosine;
        output_map[1][0] = cosine;
        output_map[1][1] = sine;
        output_map[2][2] = 1.0;
        output_map[3][3] = 1.0;

        double y0[4] = {};
        std::vector<double> predicted(static_cast<std::size_t>(4 * nvar), 0.0);
        for (int row = 0; row < 4; ++row) {
            double sum = 0.0;
            for (int i = 0; i < 6; ++i) {
                sum += output_map[row][i] * offset[static_cast<std::size_t>(i)];
            }
            y0[row] = sum;
            for (int col = 0; col < nvar; ++col) {
                double entry = 0.0;
                for (int i = 0; i < 6; ++i) {
                    entry += output_map[row][i] * map(i, col);
                }
                predicted[static_cast<std::size_t>(row * nvar + col)] = entry;
            }
        }

        for (int row = 0; row < 4; ++row) {
            const double weight = scale * stage_weight[row];
            if (weight == 0.0) {
                continue;
            }
            for (int c1 = 0; c1 < nvar; ++c1) {
                const double m1 = predicted[static_cast<std::size_t>(row * nvar + c1)];
                gradient[static_cast<std::size_t>(c1)] += 2.0 * weight * m1 * y0[row];
                for (int c2 = c1; c2 < nvar; ++c2) {
                    const double m2 = predicted[static_cast<std::size_t>(row * nvar + c2)];
                    hessian(c1, c2) += 2.0 * weight * m1 * m2;
                }
            }
        }
    };

    for (int k = 0; k < horizon; ++k) {
        Matrix Ad;
        Matrix Bd;
        const VehicleState& ref = reference[static_cast<std::size_t>(k)];
        const Input& uref = reference_input[static_cast<std::size_t>(k)];
        discrete_jacobians(ref, uref, vehicle_, dt, Ad, Bd);
        const VehicleState predicted = integrate_rk4(ref, uref, vehicle_, dt);
        const std::array<double, 6> mismatch = state_error(predicted, reference[static_cast<std::size_t>(k + 1)]);
        const std::vector<double> input_ref{uref.delta, uref.accel};
        const std::vector<double> input_effect = Bd.multiply(input_ref);

        Matrix next_sensitivity = Ad * sensitivity;
        next_sensitivity(0, 2 * k) += Bd(0, 0);
        next_sensitivity(1, 2 * k) += Bd(1, 0);
        next_sensitivity(2, 2 * k) += Bd(2, 0);
        next_sensitivity(3, 2 * k) += Bd(3, 0);
        next_sensitivity(4, 2 * k) += Bd(4, 0);
        next_sensitivity(5, 2 * k) += Bd(5, 0);
        next_sensitivity(0, 2 * k + 1) += Bd(0, 1);
        next_sensitivity(1, 2 * k + 1) += Bd(1, 1);
        next_sensitivity(2, 2 * k + 1) += Bd(2, 1);
        next_sensitivity(3, 2 * k + 1) += Bd(3, 1);
        next_sensitivity(4, 2 * k + 1) += Bd(4, 1);
        next_sensitivity(5, 2 * k + 1) += Bd(5, 1);

        std::vector<double> next_affine = Ad.multiply(affine);
        for (int i = 0; i < 6; ++i) {
            next_affine[static_cast<std::size_t>(i)] +=
                mismatch[static_cast<std::size_t>(i)] - input_effect[static_cast<std::size_t>(i)];
        }
        sensitivity = std::move(next_sensitivity);
        affine = std::move(next_affine);
        const double scale = (k + 1 == horizon) ? tuning_.terminal_scale : 1.0;
        accumulate(sensitivity, affine, reference[static_cast<std::size_t>(k + 1)].psi, scale);
    }

    const double input_weight[2] = {tuning_.r_delta, tuning_.r_accel};
    const double rate_weight[2] = {tuning_.rate_delta, tuning_.rate_accel};
    for (int k = 0; k < horizon; ++k) {
        const Input& uref = reference_input[static_cast<std::size_t>(k)];
        const double ref_u[2] = {uref.delta, uref.accel};
        for (int j = 0; j < 2; ++j) {
            const int index = 2 * k + j;
            hessian(index, index) += 2.0 * input_weight[j];
            gradient[static_cast<std::size_t>(index)] += -2.0 * input_weight[j] * ref_u[j];
        }
    }
    const double previous_u[2] = {previous_.delta, previous_.accel};
    for (int j = 0; j < 2; ++j) {
        hessian(j, j) += 2.0 * rate_weight[j];
        gradient[static_cast<std::size_t>(j)] += -2.0 * rate_weight[j] * previous_u[j];
    }
    for (int k = 1; k < horizon; ++k) {
        for (int j = 0; j < 2; ++j) {
            const int previous_index = 2 * (k - 1) + j;
            const int index = 2 * k + j;
            hessian(previous_index, previous_index) += 2.0 * rate_weight[j];
            hessian(index, index) += 2.0 * rate_weight[j];
            hessian(previous_index, index) += -2.0 * rate_weight[j];
        }
    }
    hessian.reflect_upper();

    const int nrows = 8 * horizon;
    Matrix constraints(nrows, nvar, 0.0);
    std::vector<double> upper(static_cast<std::size_t>(nrows), 0.0);
    int row = 0;
    for (int k = 0; k < horizon; ++k) {
        const int delta_index = 2 * k;
        const int accel_index = 2 * k + 1;
        constraints(row, delta_index) = 1.0;
        upper[static_cast<std::size_t>(row)] = vehicle_.delta_max;
        ++row;
        constraints(row, delta_index) = -1.0;
        upper[static_cast<std::size_t>(row)] = vehicle_.delta_max;
        ++row;
        constraints(row, accel_index) = 1.0;
        upper[static_cast<std::size_t>(row)] = vehicle_.accel_max;
        ++row;
        constraints(row, accel_index) = -1.0;
        upper[static_cast<std::size_t>(row)] = -vehicle_.accel_min;
        ++row;

        const double delta_step = vehicle_.delta_rate_max * dt;
        const double accel_step = vehicle_.accel_rate_max * dt;
        if (k == 0) {
            constraints(row, delta_index) = 1.0;
            upper[static_cast<std::size_t>(row)] = previous_.delta + delta_step;
            ++row;
            constraints(row, delta_index) = -1.0;
            upper[static_cast<std::size_t>(row)] = delta_step - previous_.delta;
            ++row;
            constraints(row, accel_index) = 1.0;
            upper[static_cast<std::size_t>(row)] = previous_.accel + accel_step;
            ++row;
            constraints(row, accel_index) = -1.0;
            upper[static_cast<std::size_t>(row)] = accel_step - previous_.accel;
            ++row;
        } else {
            const int prev_delta = 2 * (k - 1);
            const int prev_accel = prev_delta + 1;
            constraints(row, delta_index) = 1.0;
            constraints(row, prev_delta) = -1.0;
            upper[static_cast<std::size_t>(row)] = delta_step;
            ++row;
            constraints(row, prev_delta) = 1.0;
            constraints(row, delta_index) = -1.0;
            upper[static_cast<std::size_t>(row)] = delta_step;
            ++row;
            constraints(row, accel_index) = 1.0;
            constraints(row, prev_accel) = -1.0;
            upper[static_cast<std::size_t>(row)] = accel_step;
            ++row;
            constraints(row, prev_accel) = 1.0;
            constraints(row, accel_index) = -1.0;
            upper[static_cast<std::size_t>(row)] = accel_step;
            ++row;
        }
    }

    QpProblem problem;
    problem.P = std::move(hessian);
    problem.q = std::move(gradient);
    problem.G = std::move(constraints);
    problem.h = std::move(upper);

    QpSettings settings;
    settings.max_iterations = tuning_.max_iterations;
    settings.tolerance = tuning_.tolerance;

    std::vector<double> warm;
    const std::vector<double>* warm_ptr = nullptr;
    if (have_warm_ && static_cast<int>(last_solution_.size()) == nvar) {
        warm.assign(static_cast<std::size_t>(nvar), 0.0);
        for (int k = 0; k < horizon - 1; ++k) {
            warm[static_cast<std::size_t>(2 * k)] = last_solution_[static_cast<std::size_t>(2 * (k + 1))];
            warm[static_cast<std::size_t>(2 * k + 1)] = last_solution_[static_cast<std::size_t>(2 * (k + 1) + 1)];
        }
        warm[static_cast<std::size_t>(nvar - 2)] = last_solution_[static_cast<std::size_t>(nvar - 2)];
        warm[static_cast<std::size_t>(nvar - 1)] = last_solution_[static_cast<std::size_t>(nvar - 1)];
        warm_ptr = &warm;
    }

    QpSolution solution = solve_qp(problem, settings, warm_ptr);
    if (!solution.success) {
        settings.max_iterations = tuning_.max_iterations * 2;
        solution = solve_qp(problem, settings, nullptr);
    }

    Command command;
    command.iterations = solution.iterations;
    if (solution.success && static_cast<int>(solution.x.size()) == nvar) {
        command.input.delta = solution.x[0];
        command.input.accel = solution.x[1];
        command.solver_ok = true;
        last_solution_ = solution.x;
        have_warm_ = true;
    } else {
        command.input = rate_limited(reference_input[0], previous_, vehicle_, dt);
        command.solver_ok = false;
        have_warm_ = false;
    }
    previous_ = command.input;

    const auto t1 = std::chrono::steady_clock::now();
    command.solve_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    return command;
}

}  // namespace bpt
