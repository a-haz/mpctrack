#include "bpt/control.hpp"
#include "bpt/path.hpp"
#include "bpt/simulator.hpp"
#include "bpt/types.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace {

// Replays a fixed list of commands and records what it was shown, so the
// simulator bookkeeping can be checked without a real controller.
class ScriptedController final : public bpt::Controller {
public:
    explicit ScriptedController(std::vector<bpt::Command> script) : script_(std::move(script)) {}

    void reset(const bpt::VehicleState&) override { next_ = 0; }

    bpt::Command step(const bpt::VehicleState& measured) override {
        seen.push_back(measured);
        const std::size_t i = next_ < script_.size() ? next_ : script_.size() - 1;
        ++next_;
        return script_[i];
    }

    const char* name() const override { return "scripted"; }

    std::vector<bpt::VehicleState> seen;

private:
    std::vector<bpt::Command> script_;
    std::size_t next_{0};
};

bpt::Command hold(double delta, double solve_ms = 0.0) {
    bpt::Command command;
    command.input.delta = delta;
    command.input.accel = 0.0;
    command.solve_ms = solve_ms;
    return command;
}

}  // namespace

TEST(Simulator, StateOnPathUsesCurvatureSteering) {
    const double radius = 20.0;
    const bpt::Path path = bpt::make_circle(radius, bpt::kPi / 2.0);
    const bpt::VehicleParams vehicle;
    const bpt::PathPoint point = path.sample(10.0);
    const bpt::VehicleState state = bpt::state_on_path(path, 10.0, vehicle);
    EXPECT_NEAR(state.x, point.x, 1e-12);
    EXPECT_NEAR(state.y, point.y, 1e-12);
    EXPECT_NEAR(state.psi, point.psi, 1e-12);
    EXPECT_NEAR(state.delta, std::atan(vehicle.wheelbase / radius), 1e-9);
    EXPECT_EQ(state.v, 0.0);
    EXPECT_EQ(state.accel, 0.0);
}

TEST(Simulator, StraightRunHasZeroErrorAndExpectedSampleCount) {
    const bpt::Path path = bpt::make_straight_line(100.0);
    const bpt::VehicleParams vehicle;
    bpt::SimConfig config;
    config.duration = 2.0;
    config.speed0 = 8.0;
    ScriptedController controller({hold(0.0)});
    const bpt::SimResult result = bpt::simulate(controller, path, vehicle, config);

    ASSERT_EQ(result.samples.size(), 40u);
    EXPECT_NEAR(result.samples[1].t, config.dt, 1e-12);
    EXPECT_NEAR(result.samples[0].x, config.s0, 1e-12);
    EXPECT_NEAR(result.samples[1].x - result.samples[0].x, config.speed0 * config.dt, 1e-9);
    EXPECT_NEAR(result.rms_cte, 0.0, 1e-9);
    EXPECT_NEAR(result.max_heading_deg, 0.0, 1e-9);
    EXPECT_EQ(result.violations, 0);

    // With noise disabled the controller sees the true state.
    ASSERT_EQ(controller.seen.size(), result.samples.size());
    EXPECT_EQ(controller.seen[5].x, result.samples[5].x);
    EXPECT_EQ(controller.seen[5].y, result.samples[5].y);
}

TEST(Simulator, CountsSteeringRateViolationOnce) {
    const bpt::Path path = bpt::make_straight_line(100.0);
    const bpt::VehicleParams vehicle;
    bpt::SimConfig config;
    config.duration = 0.5;
    // The rate limit allows 0.4 rad/s * 0.05 s = 0.02 rad per step, so the jump to
    // 0.1 rad on step 3 is a violation and holding it afterwards is not.
    ScriptedController controller({hold(0.0), hold(0.0), hold(0.0), hold(0.1)});
    const bpt::SimResult result = bpt::simulate(controller, path, vehicle, config);

    ASSERT_EQ(result.samples.size(), 10u);
    EXPECT_EQ(result.violations, 1);
    EXPECT_FALSE(result.samples[2].violation);
    EXPECT_TRUE(result.samples[3].violation);
    EXPECT_FALSE(result.samples[4].violation);
}

TEST(Simulator, ReportsMeanAndMaxSolveTime) {
    const bpt::Path path = bpt::make_straight_line(100.0);
    bpt::SimConfig config;
    config.duration = 0.15;
    ScriptedController controller({hold(0.0, 1.0), hold(0.0, 2.0), hold(0.0, 6.0)});
    const bpt::SimResult result = bpt::simulate(controller, path, bpt::VehicleParams{}, config);

    ASSERT_EQ(result.samples.size(), 3u);
    EXPECT_NEAR(result.mean_solve_ms, 3.0, 1e-12);
    EXPECT_NEAR(result.max_solve_ms, 6.0, 1e-12);
}

TEST(Simulator, NoiseIsRepeatableForASeedAndOnlyAffectsMeasurement) {
    const bpt::Path path = bpt::make_straight_line(100.0);
    bpt::SimConfig config;
    config.duration = 1.0;
    config.noise.enabled = true;
    config.noise.seed = 7;

    ScriptedController first({hold(0.0)});
    ScriptedController second({hold(0.0)});
    const bpt::SimResult a = bpt::simulate(first, path, bpt::VehicleParams{}, config);
    const bpt::SimResult b = bpt::simulate(second, path, bpt::VehicleParams{}, config);

    ASSERT_EQ(first.seen.size(), second.seen.size());
    for (std::size_t k = 0; k < first.seen.size(); ++k) {
        EXPECT_EQ(first.seen[k].x, second.seen[k].x);
        EXPECT_EQ(first.seen[k].y, second.seen[k].y);
    }
    // The commands ignore the measurement, so the truth stays on the line.
    EXPECT_NEAR(a.rms_cte, 0.0, 1e-9);
    EXPECT_NE(first.seen[3].y, a.samples[3].y);
    EXPECT_EQ(a.rms_cte, b.rms_cte);
}
