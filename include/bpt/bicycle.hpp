#pragma once

#include "bpt/matrix.hpp"
#include "bpt/types.hpp"

namespace bpt {

// Continuous-time derivative and Jacobians of the augmented kinematic bicycle,
// including first-order steering and acceleration actuators.
std::array<double, 6> continuous_derivative(const VehicleState& state, const Input& input,
                                            const VehicleParams& params);

void continuous_jacobians(const VehicleState& state, const Input& input, const VehicleParams& params, Matrix& A,
                          Matrix& B);

// One RK4 step of the same model. The command is held constant over dt.
VehicleState integrate_rk4(const VehicleState& state, const Input& input, const VehicleParams& params,
                           double dt);

// State-transition and input matrices of that RK4 map: x+ = A x_dev + B u_dev + ...
void discrete_jacobians(const VehicleState& state, const Input& input, const VehicleParams& params,
                        double dt, Matrix& A, Matrix& B);

}  // namespace bpt
