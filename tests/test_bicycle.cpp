#include "bpt/bicycle.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace {

bpt::VehicleState shift_state(bpt::VehicleState state, int index, double delta) {
    auto values = state.to_array();
    values[static_cast<std::size_t>(index)] += delta;
    return bpt::VehicleState::from_array(values);
}

bpt::Input shift_input(bpt::Input input, int index, double delta) {
    if (index == 0) {
        input.delta += delta;
    } else {
        input.accel += delta;
    }
    return input;
}

double max_abs_diff(const bpt::Matrix& analytic, const bpt::Matrix& numeric) {
    double worst = 0.0;
    for (int r = 0; r < analytic.rows(); ++r) {
        for (int c = 0; c < analytic.cols(); ++c) {
            worst = std::max(worst, std::abs(analytic(r, c) - numeric(r, c)));
        }
    }
    return worst;
}

void expect_discrete_jacobian(const bpt::VehicleState& state, const bpt::Input& input, double dt) {
    const bpt::VehicleParams params;
    bpt::Matrix analytic_A;
    bpt::Matrix analytic_B;
    bpt::discrete_jacobians(state, input, params, dt, analytic_A, analytic_B);

    const double eps = 1e-6;
    bpt::Matrix numeric_A(6, 6, 0.0);
    bpt::Matrix numeric_B(6, 2, 0.0);
    for (int j = 0; j < 6; ++j) {
        const auto plus = bpt::integrate_rk4(shift_state(state, j, eps), input, params, dt).to_array();
        const auto minus = bpt::integrate_rk4(shift_state(state, j, -eps), input, params, dt).to_array();
        for (int i = 0; i < 6; ++i) {
            numeric_A(i, j) = (plus[static_cast<std::size_t>(i)] - minus[static_cast<std::size_t>(i)]) / (2.0 * eps);
        }
    }
    for (int j = 0; j < 2; ++j) {
        const auto plus = bpt::integrate_rk4(state, shift_input(input, j, eps), params, dt).to_array();
        const auto minus = bpt::integrate_rk4(state, shift_input(input, j, -eps), params, dt).to_array();
        for (int i = 0; i < 6; ++i) {
            numeric_B(i, j) = (plus[static_cast<std::size_t>(i)] - minus[static_cast<std::size_t>(i)]) / (2.0 * eps);
        }
    }
    EXPECT_LT(max_abs_diff(analytic_A, numeric_A), 1e-5);
    EXPECT_LT(max_abs_diff(analytic_B, numeric_B), 1e-5);
}

}  // namespace

TEST(BicycleModel, StraightConstantSpeedIsExact) {
    bpt::VehicleParams params;
    bpt::VehicleState state;
    state.v = 8.0;
    state.psi = 0.3;
    bpt::Input input;
    input.delta = state.delta;
    input.accel = 0.0;
    const double dt = 0.05;
    for (int i = 0; i < 20; ++i) {
        state = bpt::integrate_rk4(state, input, params, dt);
    }
    EXPECT_NEAR(state.x, 8.0 * std::cos(0.3), 1e-9);
    EXPECT_NEAR(state.y, 8.0 * std::sin(0.3), 1e-9);
    EXPECT_NEAR(state.psi, 0.3, 1e-12);
    EXPECT_NEAR(state.v, 8.0, 1e-12);
}

TEST(BicycleModel, StraightAccelerationIsExactWhenActuatorIsSettled) {
    bpt::VehicleParams params;
    bpt::VehicleState state;
    state.v = 5.0;
    state.accel = 1.0;
    bpt::Input input;
    input.accel = 1.0;
    const double T = 1.5;
    const double dt = 0.01;
    for (double t = 0.0; t < T - 1e-12; t += dt) {
        state = bpt::integrate_rk4(state, input, params, dt);
    }
    EXPECT_NEAR(state.v, 5.0 + T, 1e-8);
    EXPECT_NEAR(state.x, 5.0 * T + 0.5 * T * T, 1e-6);
    EXPECT_NEAR(state.y, 0.0, 1e-10);
    EXPECT_NEAR(state.psi, 0.0, 1e-12);
}

TEST(BicycleModel, ActuatorLagMatchesExponential) {
    bpt::VehicleParams params;
    bpt::VehicleState state;
    state.v = 3.0;
    bpt::Input input;
    input.delta = 0.2;
    input.accel = -1.0;
    const double T = 0.5;
    const double dt = 0.005;
    for (double t = 0.0; t < T - 1e-12; t += dt) {
        state = bpt::integrate_rk4(state, input, params, dt);
    }
    const double delta_exact = input.delta * (1.0 - std::exp(-T / params.tau_delta));
    const double accel_exact = input.accel * (1.0 - std::exp(-T / params.tau_accel));
    EXPECT_NEAR(state.delta, delta_exact, 1e-6);
    EXPECT_NEAR(state.accel, accel_exact, 1e-6);
}

TEST(BicycleModel, CircularArcMatchesConstantCurvature) {
    bpt::VehicleParams params;
    bpt::VehicleState state;
    state.v = 8.0;
    state.delta = 0.15;
    bpt::Input input;
    input.delta = 0.15;
    const double yaw_rate = state.v * std::tan(state.delta) / params.wheelbase;
    const double radius = state.v / yaw_rate;
    const double T = 2.0;
    const double dt = 0.005;
    for (double t = 0.0; t < T - 1e-12; t += dt) {
        state = bpt::integrate_rk4(state, input, params, dt);
    }
    const double psi = yaw_rate * T;
    EXPECT_NEAR(state.psi, psi, 1e-6);
    EXPECT_NEAR(state.x, radius * std::sin(psi), 1e-4);
    EXPECT_NEAR(state.y, radius * (1.0 - std::cos(psi)), 1e-4);
    EXPECT_NEAR(state.v, 8.0, 1e-8);
    EXPECT_NEAR(state.delta, 0.15, 1e-8);
}

TEST(BicycleModel, ContinuousJacobianMatchesFiniteDifference) {
    bpt::VehicleParams params;
    bpt::VehicleState state;
    state.x = 2.0;
    state.y = -1.0;
    state.psi = -0.8;
    state.v = 9.0;
    state.delta = 0.18;
    state.accel = -0.4;
    bpt::Input input{0.22, -0.7};

    bpt::Matrix analytic_A;
    bpt::Matrix analytic_B;
    bpt::continuous_jacobians(state, input, params, analytic_A, analytic_B);

    const double eps = 1e-6;
    bpt::Matrix numeric_A(6, 6, 0.0);
    bpt::Matrix numeric_B(6, 2, 0.0);
    for (int j = 0; j < 6; ++j) {
        const auto plus = bpt::continuous_derivative(shift_state(state, j, eps), input, params);
        const auto minus = bpt::continuous_derivative(shift_state(state, j, -eps), input, params);
        for (int i = 0; i < 6; ++i) {
            numeric_A(i, j) =
                (plus[static_cast<std::size_t>(i)] - minus[static_cast<std::size_t>(i)]) / (2.0 * eps);
        }
    }
    for (int j = 0; j < 2; ++j) {
        const auto plus = bpt::continuous_derivative(state, shift_input(input, j, eps), params);
        const auto minus = bpt::continuous_derivative(state, shift_input(input, j, -eps), params);
        for (int i = 0; i < 6; ++i) {
            numeric_B(i, j) =
                (plus[static_cast<std::size_t>(i)] - minus[static_cast<std::size_t>(i)]) / (2.0 * eps);
        }
    }
    EXPECT_LT(max_abs_diff(analytic_A, numeric_A), 1e-6);
    EXPECT_LT(max_abs_diff(analytic_B, numeric_B), 1e-6);
}

TEST(BicycleModel, DiscreteJacobianStraight) {
    bpt::VehicleState state;
    state.psi = 0.4;
    state.v = 8.0;
    state.accel = 0.2;
    bpt::Input input{0.0, 0.2};
    expect_discrete_jacobian(state, input, 0.05);
}

TEST(BicycleModel, DiscreteJacobianTurning) {
    bpt::VehicleState state;
    state.psi = -1.1;
    state.v = 11.0;
    state.delta = 0.2;
    state.accel = -0.3;
    bpt::Input input{0.25, -0.5};
    expect_discrete_jacobian(state, input, 0.05);
}

TEST(BicycleModel, DiscreteJacobianSlowTightSteer) {
    bpt::VehicleState state;
    state.psi = 2.4;
    state.v = 1.5;
    state.delta = -0.35;
    state.accel = 0.4;
    bpt::Input input{-0.3, 0.6};
    expect_discrete_jacobian(state, input, 0.02);
}

TEST(BicycleModel, DiscreteJacobianMatchesEulerForSmallStep) {
    bpt::VehicleParams params;
    bpt::VehicleState state;
    state.psi = 0.5;
    state.v = 7.0;
    state.delta = 0.1;
    bpt::Input input{0.1, -0.2};
    const double dt = 1e-4;
    bpt::Matrix Ad;
    bpt::Matrix Bd;
    bpt::Matrix Ac;
    bpt::Matrix Bc;
    bpt::discrete_jacobians(state, input, params, dt, Ad, Bd);
    bpt::continuous_jacobians(state, input, params, Ac, Bc);
    for (int r = 0; r < 6; ++r) {
        for (int c = 0; c < 6; ++c) {
            const double euler = (r == c ? 1.0 : 0.0) + dt * Ac(r, c);
            EXPECT_NEAR(Ad(r, c), euler, 1e-6);
        }
        for (int c = 0; c < 2; ++c) {
            EXPECT_NEAR(Bd(r, c), dt * Bc(r, c), 5e-6);
        }
    }
}
