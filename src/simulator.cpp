#include "bpt/simulator.hpp"

#include "bpt/bicycle.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace bpt {

VehicleState state_on_path(const Path& path, double s, const VehicleParams& params) {
    const PathPoint point = path.sample(s);
    VehicleState state;
    state.x = point.x;
    state.y = point.y;
    state.psi = point.psi;
    state.v = 0.0;
    state.delta = std::clamp(std::atan(params.wheelbase * point.kappa), -params.delta_max, params.delta_max);
    state.accel = 0.0;
    return state;
}

SimResult simulate(Controller& controller, const Path& path, const VehicleParams& params, const SimConfig& config) {
    VehicleState truth = state_on_path(path, config.s0, params);
    truth.v = config.speed0;
    controller.reset(truth);

    Input previous{truth.delta, truth.accel};
    double s_hint = -1.0;
    std::mt19937 rng(config.noise.seed);
    std::normal_distribution<double> noise_xy(0.0, config.noise.sigma_xy);
    std::normal_distribution<double> noise_psi(0.0, config.noise.sigma_psi);
    std::normal_distribution<double> noise_v(0.0, config.noise.sigma_v);

    const int steps = std::max(1, static_cast<int>(std::llround(config.duration / config.dt)));
    SimResult result;
    result.samples.reserve(static_cast<std::size_t>(steps));

    double sum_cte2 = 0.0;
    double sum_heading2 = 0.0;
    double sum_solve = 0.0;

    for (int k = 0; k < steps; ++k) {
        VehicleState measured = truth;
        if (config.noise.enabled) {
            measured.x += noise_xy(rng);
            measured.y += noise_xy(rng);
            measured.psi = wrap_to_pi(measured.psi + noise_psi(rng));
            measured.v += noise_v(rng);
        }

        const Command command = controller.step(measured);
        const bool violation =
            command_violates(command.input, previous, params, config.dt, config.violation_tolerance);

        const Path::Projection projection = path.project(truth.x, truth.y, s_hint, 40.0);
        s_hint = projection.s;

        Sample sample;
        sample.t = static_cast<double>(k) * config.dt;
        sample.x = truth.x;
        sample.y = truth.y;
        sample.psi = truth.psi;
        sample.v = truth.v;
        sample.cte = projection.cte;
        sample.heading_deg = degrees(wrap_to_pi(truth.psi - projection.psi));
        sample.delta_cmd = command.input.delta;
        sample.accel_cmd = command.input.accel;
        sample.solve_ms = command.solve_ms;
        sample.violation = violation;
        result.samples.push_back(sample);

        sum_cte2 += sample.cte * sample.cte;
        sum_heading2 += sample.heading_deg * sample.heading_deg;
        sum_solve += sample.solve_ms;
        result.max_cte = std::max(result.max_cte, std::abs(sample.cte));
        result.max_heading_deg = std::max(result.max_heading_deg, std::abs(sample.heading_deg));
        result.max_solve_ms = std::max(result.max_solve_ms, sample.solve_ms);
        if (violation) {
            ++result.violations;
        }

        truth = integrate_rk4(truth, command.input, params, config.dt);
        previous = command.input;
    }

    const double count = static_cast<double>(result.samples.size());
    result.rms_cte = std::sqrt(sum_cte2 / count);
    result.rms_heading_deg = std::sqrt(sum_heading2 / count);
    result.mean_solve_ms = sum_solve / count;
    return result;
}

}  // namespace bpt
