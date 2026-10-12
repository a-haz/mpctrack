#include "bpt/control.hpp"
#include "bpt/simulator.hpp"

#include <gtest/gtest.h>

#include <cmath>

TEST(Constraints, MpcRespectsSteerAccelAndRateLimits) {
    bpt::CurvePathSpec spec;
    spec.straight_lead = 12.0;
    spec.transition = 6.0;
    spec.radius = 8.0;
    spec.angle = 1.0;
    spec.straight_tail = 20.0;
    spec.sample_ds = 0.05;
    const bpt::Path path = bpt::make_curve_path(spec, "tight_curve");

    const bpt::VehicleParams params;
    bpt::MpcParams tuning;
    tuning.horizon = 16;
    tuning.w_cte = 12.0;
    tuning.r_delta = 0.05;
    const double speed = 12.0;
    bpt::LtvMpc controller(path, params, tuning, speed);

    bpt::SimConfig config;
    config.dt = tuning.dt;
    config.duration = 3.0;
    config.s0 = 2.0;
    config.speed0 = speed;
    config.noise.enabled = false;
    config.violation_tolerance = 1e-4;

    const bpt::SimResult result = bpt::simulate(controller, path, params, config);
    EXPECT_EQ(result.violations, 0);
    ASSERT_FALSE(result.samples.empty());

    double max_abs_delta = 0.0;
    double max_abs_delta_step = 0.0;
    double max_abs_accel_step = 0.0;
    double previous_delta = result.samples.front().delta_cmd;
    double previous_accel = result.samples.front().accel_cmd;
    for (std::size_t i = 0; i < result.samples.size(); ++i) {
        const bpt::Sample& sample = result.samples[i];
        EXPECT_LE(sample.delta_cmd, params.delta_max + 1e-4);
        EXPECT_GE(sample.delta_cmd, -params.delta_max - 1e-4);
        EXPECT_LE(sample.accel_cmd, params.accel_max + 1e-4);
        EXPECT_GE(sample.accel_cmd, params.accel_min - 1e-4);
        max_abs_delta = std::max(max_abs_delta, std::abs(sample.delta_cmd));
        if (i > 0) {
            max_abs_delta_step = std::max(max_abs_delta_step, std::abs(sample.delta_cmd - previous_delta));
            max_abs_accel_step = std::max(max_abs_accel_step, std::abs(sample.accel_cmd - previous_accel));
        }
        previous_delta = sample.delta_cmd;
        previous_accel = sample.accel_cmd;
    }
    EXPECT_GT(max_abs_delta, 0.05);
    EXPECT_LE(max_abs_delta_step, params.delta_rate_max * config.dt + 1e-4);
    EXPECT_LE(max_abs_accel_step, params.accel_rate_max * config.dt + 1e-4);
}
