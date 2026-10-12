#include "bpt/control.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

namespace bpt {

PurePursuit::PurePursuit(const Path& path, VehicleParams vehicle, PurePursuitParams tuning, double speed_ref,
                         double dt)
    : path_(path), vehicle_(vehicle), tuning_(tuning), speed_ref_(speed_ref), dt_(dt) {}

void PurePursuit::reset(const VehicleState& state) {
    s_hint_ = -1.0;
    integral_ = 0.0;
    previous_speed_error_ = speed_ref_ - state.v;
}

Command PurePursuit::step(const VehicleState& measured) {
    const auto t0 = std::chrono::steady_clock::now();
    const Path::Projection proj = path_.project(measured.x, measured.y, s_hint_, 25.0);
    s_hint_ = proj.s;

    const double speed_for_lookahead = std::max(measured.v, 0.0);
    const double lookahead = std::clamp(tuning_.lookahead_gain * speed_for_lookahead, tuning_.lookahead_min,
                                        tuning_.lookahead_max);

    PathPoint target = path_.sample(proj.s);
    const double s_limit = path_.length();
    for (double s = proj.s; s <= s_limit; s += 0.2) {
        const PathPoint candidate = path_.sample(s);
        target = candidate;
        const double dist = std::hypot(candidate.x - measured.x, candidate.y - measured.y);
        if (dist >= lookahead) {
            break;
        }
    }

    const double alpha =
        wrap_to_pi(std::atan2(target.y - measured.y, target.x - measured.x) - measured.psi);
    double delta = std::atan2(2.0 * vehicle_.wheelbase * std::sin(alpha), std::max(lookahead, 0.5));
    delta = std::clamp(delta, -vehicle_.delta_max, vehicle_.delta_max);

    const double speed_error = speed_ref_ - measured.v;
    integral_ += speed_error * dt_;
    integral_ = std::clamp(integral_, -tuning_.integral_limit, tuning_.integral_limit);
    const double derivative = (speed_error - previous_speed_error_) / std::max(dt_, 1e-3);
    previous_speed_error_ = speed_error;
    double accel = tuning_.speed_kp * speed_error + tuning_.speed_ki * integral_ + tuning_.speed_kd * derivative;
    accel = std::clamp(accel, vehicle_.accel_min, vehicle_.accel_max);

    const auto t1 = std::chrono::steady_clock::now();
    Command command;
    command.input.delta = delta;
    command.input.accel = accel;
    command.solve_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    command.iterations = 0;
    command.solver_ok = true;
    return command;
}

}  // namespace bpt
