#include "bpt/qp_solver.hpp"

#include <gtest/gtest.h>

#include <cmath>

namespace {

bpt::QpSolution solve_or_fail(const bpt::QpProblem& problem) {
    bpt::QpSettings settings;
    settings.tolerance = 1e-9;
    settings.max_iterations = 80;
    const bpt::QpSolution solution = bpt::solve_qp(problem, settings);
    EXPECT_TRUE(solution.success);
    return solution;
}

}  // namespace

TEST(QpSolver, UnconstrainedDiagonal) {
    bpt::QpProblem problem;
    problem.P = bpt::Matrix(2, 2, 0.0);
    problem.P(0, 0) = 2.0;
    problem.P(1, 1) = 4.0;
    problem.q = {1.0, -2.0};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 2u);
    EXPECT_NEAR(solution.x[0], -0.5, 1e-8);
    EXPECT_NEAR(solution.x[1], 0.5, 1e-8);
}

TEST(QpSolver, UpperBoundActive) {
    // 1/2 x^2 - 3x subject to x <= 1. Minimizer of the parabola is 3, so x = 1.
    bpt::QpProblem problem;
    problem.P = bpt::Matrix(1, 1, 1.0);
    problem.q = {-3.0};
    problem.G = bpt::Matrix(1, 1, 1.0);
    problem.h = {1.0};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 1u);
    EXPECT_NEAR(solution.x[0], 1.0, 1e-6);
    EXPECT_LE(solution.x[0], 1.0 + 1e-8);
}

TEST(QpSolver, TwoSidedBoundSelectsUpperEdge) {
    bpt::QpProblem problem;
    problem.P = bpt::Matrix(1, 1, 1.0);
    problem.q = {-3.0};
    problem.G = bpt::Matrix(2, 1, 0.0);
    problem.G(0, 0) = 1.0;
    problem.G(1, 0) = -1.0;
    problem.h = {1.0, 1.0};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 1u);
    EXPECT_NEAR(solution.x[0], 1.0, 1e-6);
}

TEST(QpSolver, LooseBoundsLeaveInteriorMinimizer) {
    // 1/2 * 2 * (x - 0.2)^2 => P = 2, q = -0.4, bounds [-2, 2].
    bpt::QpProblem problem;
    problem.P = bpt::Matrix(1, 1, 2.0);
    problem.q = {-0.4};
    problem.G = bpt::Matrix(2, 1, 0.0);
    problem.G(0, 0) = 1.0;
    problem.G(1, 0) = -1.0;
    problem.h = {2.0, 2.0};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 1u);
    EXPECT_NEAR(solution.x[0], 0.2, 1e-6);
}

TEST(QpSolver, EqualityConstraint) {
    // P = [[4, 1], [1, 2]], q = 0, x + y = 1 -> (0.25, 0.75).
    bpt::QpProblem problem;
    problem.P = bpt::Matrix(2, 2, 0.0);
    problem.P(0, 0) = 4.0;
    problem.P(0, 1) = 1.0;
    problem.P(1, 0) = 1.0;
    problem.P(1, 1) = 2.0;
    problem.q = {0.0, 0.0};
    problem.A = bpt::Matrix(1, 2, 1.0);
    problem.b = {1.0};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 2u);
    EXPECT_NEAR(solution.x[0], 0.25, 1e-8);
    EXPECT_NEAR(solution.x[1], 0.75, 1e-8);
    EXPECT_NEAR(solution.x[0] + solution.x[1], 1.0, 1e-8);
}

TEST(QpSolver, HalfPlaneInequality) {
    // 1/2 ||x||^2 - 2 x - 2 y subject to x + y <= 1 -> (0.5, 0.5).
    bpt::QpProblem problem;
    problem.P = bpt::Matrix::identity(2);
    problem.q = {-2.0, -2.0};
    problem.G = bpt::Matrix(1, 2, 1.0);
    problem.h = {1.0};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 2u);
    EXPECT_NEAR(solution.x[0], 0.5, 1e-6);
    EXPECT_NEAR(solution.x[1], 0.5, 1e-6);
    EXPECT_LE(solution.x[0] + solution.x[1], 1.0 + 1e-7);
}

TEST(QpSolver, BoxAndRateLimits) {
    // Targets (0, 1) with |u1 - u0| <= 0.2 and loose boxes. Solution is (0.4, 0.6).
    bpt::QpProblem problem;
    problem.P = bpt::Matrix(2, 2, 0.0);
    problem.P(0, 0) = 2.0;
    problem.P(1, 1) = 2.0;
    problem.q = {0.0, -2.0};
    problem.G = bpt::Matrix(6, 2, 0.0);
    problem.h = {2.0, 2.0, 2.0, 2.0, 0.2, 0.2};
    problem.G(0, 0) = 1.0;
    problem.G(1, 0) = -1.0;
    problem.G(2, 1) = 1.0;
    problem.G(3, 1) = -1.0;
    problem.G(4, 1) = 1.0;
    problem.G(4, 0) = -1.0;
    problem.G(5, 0) = 1.0;
    problem.G(5, 1) = -1.0;
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 2u);
    EXPECT_NEAR(solution.x[0], 0.4, 1e-5);
    EXPECT_NEAR(solution.x[1], 0.6, 1e-5);
    EXPECT_LE(std::abs(solution.x[1] - solution.x[0]), 0.2 + 1e-6);
}

TEST(QpSolver, EqualityWithActiveInequality) {
    // Minimize x^2 + y^2 (P = 2 I) subject to x + y = 1 and x >= 0.8.
    bpt::QpProblem problem;
    problem.P = bpt::Matrix::identity(2);
    problem.P(0, 0) = 2.0;
    problem.P(1, 1) = 2.0;
    problem.q = {0.0, 0.0};
    problem.A = bpt::Matrix(1, 2, 1.0);
    problem.b = {1.0};
    problem.G = bpt::Matrix(1, 2, 0.0);
    problem.G(0, 0) = -1.0;
    problem.h = {-0.8};
    const bpt::QpSolution solution = solve_or_fail(problem);
    ASSERT_EQ(solution.x.size(), 2u);
    EXPECT_NEAR(solution.x[0], 0.8, 1e-5);
    EXPECT_NEAR(solution.x[1], 0.2, 1e-5);
}
