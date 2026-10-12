#include "bpt/bicycle.hpp"
#include "bpt/control.hpp"
#include "bpt/simulator.hpp"

#include <gtest/gtest.h>

#include <cmath>

namespace {

bpt::SimConfig straight_config() {
    bpt::SimConfig config;
    config.dt = 0.05;
    config.duration = 8.0;
    config.s0 = 0.0;
    config.speed0 = 8.0;
    config.noise.enabled = false;
    return config;
}

bpt::VehicleState offset_start() {
    bpt::VehicleState state;
    state.x = 0.0;
    state.y = 1.2;
    state.psi = 0.12;
    state.v = 8.0;
    return state;
}

void expect_straight_recovery(bpt::Controller& controller, const bpt::Path& path) {
    const bpt::VehicleParams params;
    bpt::SimConfig config = straight_config();
    bpt::VehicleState truth = offset_start();
    controller.reset(truth);
    bpt::Input previous{};
    double s_hint = -1.0;
    const int steps = static_cast<int>(std::llround(config.duration / config.dt));
    double initial_cte = 0.0;
    double final_cte = 0.0;
    double final_heading = 0.0;
    for (int k = 0; k < steps; ++k) {
        const bpt::Command command = controller.step(truth);
        EXPECT_TRUE(std::isfinite(command.input.delta));
        EXPECT_TRUE(std::isfinite(command.input.accel));
        const auto projection = path.project(truth.x, truth.y, s_hint, 30.0);
        s_hint = projection.s;
        if (k == 0) {
            initial_cte = std::abs(projection.cte);
        }
        final_cte = projection.cte;
        final_heading = bpt::degrees(bpt::wrap_to_pi(truth.psi - projection.psi));
        truth = bpt::integrate_rk4(truth, command.input, params, config.dt);
        previous = command.input;
        (void)previous;
    }
    EXPECT_GT(initial_cte, 1.0);
    EXPECT_LT(std::abs(final_cte), 0.20);
    EXPECT_LT(std::abs(final_heading), 5.0);
    EXPECT_LT(std::abs(final_cte), 0.5 * initial_cte);
}

}  // namespace

TEST(PathGeometry, StraightCrossTrackSign) {
    const bpt::Path path = bpt::make_straight_line(100.0);
    const auto left = path.project(10.0, 1.5, -1.0, 10.0);
    const auto right = path.project(10.0, -0.4, -1.0, 10.0);
    EXPECT_NEAR(left.cte, 1.5, 1e-9);
    EXPECT_NEAR(right.cte, -0.4, 1e-9);
    EXPECT_NEAR(left.psi, 0.0, 1e-12);
}

TEST(PathGeometry, CircleCurvature) {
    const double radius = 30.0;
    const bpt::Path path = bpt::make_circle(radius, bpt::kPi / 2.0);
    const bpt::PathPoint mid = path.sample(path.nominal_length * 0.5);
    EXPECT_NEAR(mid.kappa, 1.0 / radius, 1e-9);
    EXPECT_NEAR(mid.psi, bpt::kPi / 4.0, 1e-6);
}

TEST(PathGeometry, StraightToCurveReachesDesignRadius) {
    const bpt::Path path = bpt::make_straight_to_curve();
    bool found = false;
    for (const bpt::PathPoint& point : path.points) {
        if (point.s < path.nominal_length && std::abs(point.kappa - 1.0 / 22.0) < 1e-3) {
            found = true;
            break;
        }
    }
    EXPECT_TRUE(found);
    EXPECT_GT(path.length(), path.nominal_length);
}

TEST(PathGeometry, FigureEightIsClosedAndCurved) {
    const bpt::Path path = bpt::make_figure_eight();
    ASSERT_FALSE(path.points.empty());
    double max_kappa = 0.0;
    for (const bpt::PathPoint& point : path.points) {
        max_kappa = std::max(max_kappa, std::abs(point.kappa));
        EXPECT_TRUE(std::isfinite(point.kappa));
        EXPECT_TRUE(std::isfinite(point.psi));
    }
    EXPECT_GT(path.nominal_length, 200.0);
    // Peak curvature of this lemniscate stays inside the steering limit.
    EXPECT_GT(max_kappa, 0.04);
    EXPECT_LT(max_kappa, 0.12);
    EXPECT_GT(path.length(), 2.0 * path.nominal_length);
}

TEST(ClosedLoop, MpcTracksStraightLineFromOffset) {
    const bpt::Path path = bpt::make_straight_line(250.0);
    bpt::LtvMpc controller(path, bpt::VehicleParams{}, bpt::MpcParams{}, 8.0);
    expect_straight_recovery(controller, path);
}

TEST(ClosedLoop, PurePursuitTracksStraightLineFromOffset) {
    const bpt::Path path = bpt::make_straight_line(250.0);
    bpt::PurePursuit controller(path, bpt::VehicleParams{}, bpt::PurePursuitParams{}, 8.0, 0.05);
    expect_straight_recovery(controller, path);
}

TEST(ClosedLoop, MpcHoldsZeroCommandWhenAligned) {
    const bpt::Path path = bpt::make_straight_line(200.0);
    const bpt::VehicleParams params;
    bpt::LtvMpc controller(path, params, bpt::MpcParams{}, 10.0);
    bpt::VehicleState state = bpt::state_on_path(path, 5.0, params);
    state.v = 10.0;
    controller.reset(state);
    const bpt::Command command = controller.step(state);
    EXPECT_TRUE(command.solver_ok);
    EXPECT_NEAR(command.input.delta, 0.0, 1e-3);
    EXPECT_NEAR(command.input.accel, 0.0, 1e-3);
}

TEST(Angles, WrapToPi) {
    EXPECT_NEAR(bpt::wrap_to_pi(0.0), 0.0, 1e-15);
    EXPECT_NEAR(std::abs(bpt::wrap_to_pi(3.0 * bpt::kPi)), bpt::kPi, 1e-12);
    EXPECT_NEAR(bpt::wrap_to_pi(-3.0 * bpt::kPi / 2.0), bpt::kPi / 2.0, 1e-12);
}
