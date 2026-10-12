#pragma once

#include <array>
#include <cmath>
#include <string>

namespace bpt {

constexpr double kPi = 3.141592653589793238462643383279502884;

inline double wrap_to_pi(double angle) {
    const double two_pi = 2.0 * kPi;
    angle = std::fmod(angle + kPi, two_pi);
    if (angle < 0.0) {
        angle += two_pi;
    }
    return angle - kPi;
}

inline double degrees(double radians) { return radians * (180.0 / kPi); }

inline double radians(double degrees_value) { return degrees_value * (kPi / 180.0); }

struct VehicleParams {
    double wheelbase{2.7};
    double delta_max{0.50};
    double accel_max{2.0};
    double accel_min{-3.5};
    double delta_rate_max{0.40};
    double accel_rate_max{5.0};
    double tau_delta{0.15};
    double tau_accel{0.18};
};

struct VehicleState {
    double x{0.0};
    double y{0.0};
    double psi{0.0};
    double v{0.0};
    double delta{0.0};
    double accel{0.0};

    std::array<double, 6> to_array() const { return {x, y, psi, v, delta, accel}; }

    static VehicleState from_array(const std::array<double, 6>& a) {
        VehicleState s;
        s.x = a[0];
        s.y = a[1];
        s.psi = a[2];
        s.v = a[3];
        s.delta = a[4];
        s.accel = a[5];
        return s;
    }
};

struct Input {
    double delta{0.0};
    double accel{0.0};

    std::array<double, 2> to_array() const { return {delta, accel}; }
};

inline bool command_violates(const Input& u, const Input& previous, const VehicleParams& params,
                             double dt, double tolerance) {
    if (u.delta > params.delta_max + tolerance || u.delta < -params.delta_max - tolerance) {
        return true;
    }
    if (u.accel > params.accel_max + tolerance || u.accel < params.accel_min - tolerance) {
        return true;
    }
    const double delta_step = params.delta_rate_max * dt + tolerance;
    const double accel_step = params.accel_rate_max * dt + tolerance;
    if (std::abs(u.delta - previous.delta) > delta_step) {
        return true;
    }
    if (std::abs(u.accel - previous.accel) > accel_step) {
        return true;
    }
    return false;
}

}  // namespace bpt
