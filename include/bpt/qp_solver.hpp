#pragma once

#include "bpt/matrix.hpp"

#include <vector>

namespace bpt {

// Dense convex quadratic program
//
//   minimize    1/2 x^T P x + q^T x
//   subject to  G x <= h
//               A x  = b
//
// P must be symmetric positive definite when inequalities are present.
// Equality-only and unconstrained problems are solved by a direct KKT factorization.
// Inequalities use a Mehrotra predictor-corrector interior-point method.

struct QpProblem {
    Matrix P;
    std::vector<double> q;
    Matrix G;
    std::vector<double> h;
    Matrix A;
    std::vector<double> b;
};

struct QpSettings {
    int max_iterations{60};
    double tolerance{1e-8};
    double fraction_to_boundary{0.995};
};

struct QpSolution {
    std::vector<double> x;
    int iterations{0};
    bool success{false};
    double primal_residual{0.0};
    double dual_residual{0.0};
    double complementarity{0.0};
};

QpSolution solve_qp(const QpProblem& problem, const QpSettings& settings = {},
                    const std::vector<double>* warm_start = nullptr);

}  // namespace bpt
