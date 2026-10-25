#include "bpt/bicycle.hpp"
#include "bpt/control.hpp"
#include "bpt/path.hpp"
#include "bpt/types.hpp"

#include <gtest/gtest.h>

#include <cmath>

namespace {

bpt::VehicleState straight_state(double y, double psi, double v) {
    bpt::VehicleState state;
    state.x = 10.0;
    state.y = y;
    state.psi = psi;
    state.v = v;
    return state;
}

}  // namespace

TEST(PurePursuit, ZeroSteerOnStraightPathWithoutOffset) {
    const bpt::Path path = bpt::make_straight_line(200.0);
    bpt::PurePursuit controller(path, bpt::VehicleParams{}, bpt::PurePursuitParams{}, 8.0, 0.05);
    const bpt::VehicleState state = straight_state(0.0, 0.0, 8.0);
    controller.reset(state);
    const bpt::Command command = controller.step(state);
    EXPECT_NEAR(command.input.delta, 0.0, 1e-9);
    EXPECT_NEAR(command.input.accel, 0.0, 1e-9);
    EXPECT_TRUE(command.solver_ok);
}

// Positive delta gives a positive yaw rate (psi_dot = v tan(delta) / L), so it turns left.
TEST(PurePursuit, OffsetRightSteersLeft) {
    const bpt::Path path = bpt::make_straight_line(200.0);
    const bpt::VehicleParams vehicle;
    bpt::PurePursuit controller(path, vehicle, bpt::PurePursuitParams{}, 8.0, 0.05);
    const bpt::VehicleState right = straight_state(-1.0, 0.0, 8.0);
    controller.reset(right);
    const bpt::Command command = controller.step(right);
    EXPECT_GT(command.input.delta, 0.01);
    EXPECT_NEAR(path.project(right.x, right.y, -1.0, 10.0).cte, -1.0, 1e-9);

    const bpt::VehicleState stepped = bpt::integrate_rk4(right, command.input, vehicle, 0.5);
    EXPECT_GT(stepped.psi, 0.0);
}

TEST(PurePursuit, OffsetLeftSteersRight) {
    const bpt::Path path = bpt::make_straight_line(200.0);
    bpt::PurePursuit controller(path, bpt::VehicleParams{}, bpt::PurePursuitParams{}, 8.0, 0.05);
    const bpt::VehicleState left = straight_state(1.0, 0.0, 8.0);
    controller.reset(left);
    EXPECT_LT(controller.step(left).input.delta, -0.01);
}

TEST(PurePursuit, CommandsStayWithinVehicleLimits) {
    const bpt::Path path = bpt::make_straight_line(200.0);
    const bpt::VehicleParams vehicle;
    const double offsets[] = {-30.0, -3.0, 3.0, 30.0};
    const double headings[] = {-1.5, 0.0, 1.5};
    const double speeds[] = {0.0, 8.0, 30.0};
    for (double y : offsets) {
        for (double psi : headings) {
            for (double v : speeds) {
                bpt::PurePursuit controller(path, vehicle, bpt::PurePursuitParams{}, 10.0, 0.05);
                const bpt::VehicleState state = straight_state(y, psi, v);
                controller.reset(state);
                const bpt::Command command = controller.step(state);
                EXPECT_LE(std::abs(command.input.delta), vehicle.delta_max + 1e-12);
                EXPECT_LE(command.input.accel, vehicle.accel_max + 1e-12);
                EXPECT_GE(command.input.accel, vehicle.accel_min - 1e-12);
            }
        }
    }

    bpt::PurePursuit far(path, vehicle, bpt::PurePursuitParams{}, 10.0, 0.05);
    const bpt::VehicleState far_right = straight_state(-30.0, 0.0, 0.0);
    far.reset(far_right);
    const bpt::Command saturated = far.step(far_right);
    EXPECT_NEAR(saturated.input.delta, vehicle.delta_max, 1e-12);
    EXPECT_NEAR(saturated.input.accel, vehicle.accel_max, 1e-12);
}

TEST(PurePursuit, ClosedLoopOnFigureEightStaysWithinLimits) {
    const bpt::Path path = bpt::make_figure_eight();
    const bpt::VehicleParams vehicle;
    const double dt = 0.05;
    bpt::PurePursuit controller(path, vehicle, bpt::PurePursuitParams{}, 8.0, dt);
    bpt::VehicleState state;
    state.x = path.points.front().x;
    state.y = path.points.front().y;
    state.psi = path.points.front().psi;
    state.v = 8.0;
    controller.reset(state);
    for (int k = 0; k < 400; ++k) {
        const bpt::Command command = controller.step(state);
        ASSERT_TRUE(std::isfinite(command.input.delta));
        ASSERT_TRUE(std::isfinite(command.input.accel));
        EXPECT_LE(std::abs(command.input.delta), vehicle.delta_max + 1e-12);
        EXPECT_LE(command.input.accel, vehicle.accel_max + 1e-12);
        EXPECT_GE(command.input.accel, vehicle.accel_min - 1e-12);
        state = bpt::integrate_rk4(state, command.input, vehicle, dt);
    }
}

TEST(PurePursuit, SpeedLoopReachesReference) {
    const bpt::Path path = bpt::make_straight_line(400.0);
    const bpt::VehicleParams vehicle;
    const double dt = 0.05;
    const double speed_ref = 10.0;
    bpt::PurePursuit controller(path, vehicle, bpt::PurePursuitParams{}, speed_ref, dt);
    bpt::VehicleState state = straight_state(0.0, 0.0, 3.0);
    state.x = 0.0;
    controller.reset(state);
    for (int k = 0; k < 400; ++k) {
        const bpt::Command command = controller.step(state);
        state = bpt::integrate_rk4(state, command.input, vehicle, dt);
    }
    EXPECT_NEAR(state.v, speed_ref, 0.05);
    EXPECT_NEAR(state.accel, 0.0, 0.05);
    EXPECT_NEAR(state.y, 0.0, 1e-6);
}

TEST(PurePursuit, SpeedLoopSlowsDownToReference) {
    const bpt::Path path = bpt::make_straight_line(400.0);
    const bpt::VehicleParams vehicle;
    const double dt = 0.05;
    const double speed_ref = 6.0;
    bpt::PurePursuit controller(path, vehicle, bpt::PurePursuitParams{}, speed_ref, dt);
    bpt::VehicleState state = straight_state(0.0, 0.0, 14.0);
    state.x = 0.0;
    controller.reset(state);
    for (int k = 0; k < 400; ++k) {
        const bpt::Command command = controller.step(state);
        state = bpt::integrate_rk4(state, command.input, vehicle, dt);
    }
    EXPECT_NEAR(state.v, speed_ref, 0.05);
}
