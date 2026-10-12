#pragma once

#include "bpt/control.hpp"

#include <vector>

namespace bpt {

struct NoiseConfig {
    bool enabled{false};
    double sigma_xy{0.025};
    double sigma_psi{0.004};
    double sigma_v{0.05};
    unsigned seed{1};
};

struct SimConfig {
    double dt{0.05};
    double duration{10.0};
    double s0{8.0};
    double speed0{8.0};
    double violation_tolerance{1e-4};
    NoiseConfig noise{};
};

struct Sample {
    double t{0.0};
    double x{0.0};
    double y{0.0};
    double psi{0.0};
    double v{0.0};
    double cte{0.0};
    double heading_deg{0.0};
    double delta_cmd{0.0};
    double accel_cmd{0.0};
    double solve_ms{0.0};
    bool violation{false};
};

struct SimResult {
    std::vector<Sample> samples;
    int violations{0};
    double rms_cte{0.0};
    double max_cte{0.0};
    double rms_heading_deg{0.0};
    double max_heading_deg{0.0};
    double mean_solve_ms{0.0};
    double max_solve_ms{0.0};
};

VehicleState state_on_path(const Path& path, double s, const VehicleParams& params);

SimResult simulate(Controller& controller, const Path& path, const VehicleParams& params,
                   const SimConfig& config);

}  // namespace bpt
