#pragma once

#include "bpt/path.hpp"
#include "bpt/types.hpp"

#include <vector>

namespace bpt {

struct Command {
    Input input;
    double solve_ms{0.0};
    int iterations{0};
    bool solver_ok{true};
};

struct MpcParams {
    int horizon{20};
    double dt{0.05};
    double w_cte{8.0};
    double w_along{0.5};
    double w_heading{6.0};
    double w_speed{1.2};
    double r_delta{0.25};
    double r_accel{0.10};
    double rate_delta{0.80};
    double rate_accel{0.08};
    double terminal_scale{5.0};
    int max_iterations{35};
    double tolerance{1e-6};
};

struct PurePursuitParams {
    double lookahead_gain{0.70};
    double lookahead_min{2.5};
    double lookahead_max{12.0};
    double speed_kp{1.5};
    double speed_ki{0.20};
    double speed_kd{0.0};
    double integral_limit{1.5};
};

class Controller {
public:
    virtual ~Controller() = default;
    virtual void reset(const VehicleState& state) = 0;
    virtual Command step(const VehicleState& measured) = 0;
    virtual const char* name() const = 0;
};

class LtvMpc final : public Controller {
public:
    LtvMpc(const Path& path, VehicleParams vehicle, MpcParams tuning, double speed_ref);

    void reset(const VehicleState& state) override;
    Command step(const VehicleState& measured) override;
    const char* name() const override { return "ltv_mpc"; }

    const std::vector<double>& last_solution() const { return last_solution_; }

private:
    const Path& path_;
    VehicleParams vehicle_;
    MpcParams tuning_;
    double speed_ref_;
    double s_hint_{-1.0};
    Input previous_{};
    std::vector<double> last_solution_;
    bool have_warm_{false};
};

class PurePursuit final : public Controller {
public:
    PurePursuit(const Path& path, VehicleParams vehicle, PurePursuitParams tuning, double speed_ref,
                double dt);

    void reset(const VehicleState& state) override;
    Command step(const VehicleState& measured) override;
    const char* name() const override { return "pure_pursuit"; }

private:
    const Path& path_;
    VehicleParams vehicle_;
    PurePursuitParams tuning_;
    double speed_ref_;
    double dt_;
    double s_hint_{-1.0};
    double integral_{0.0};
    double previous_speed_error_{0.0};
};

}  // namespace bpt
