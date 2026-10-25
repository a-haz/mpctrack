#include "bpt/path.hpp"
#include "bpt/types.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>

namespace {

bpt::Path two_point_path() {
    bpt::Path path;
    path.name = "two_point";
    path.points.push_back(bpt::PathPoint{0.0, 0.0, 0.0, 0.0, 0.0});
    path.points.push_back(bpt::PathPoint{2.0, 2.0, 1.0, 0.4, 0.2});
    path.nominal_length = path.length();
    return path;
}

}  // namespace

TEST(PathSample, InterpolatesBetweenPoints) {
    const bpt::Path path = two_point_path();
    const bpt::PathPoint p = path.sample(0.5);
    EXPECT_NEAR(p.s, 0.5, 1e-12);
    EXPECT_NEAR(p.x, 0.5, 1e-12);
    EXPECT_NEAR(p.y, 0.25, 1e-12);
    EXPECT_NEAR(p.psi, 0.1, 1e-12);
    EXPECT_NEAR(p.kappa, 0.05, 1e-12);

    const bpt::PathPoint mid = path.sample(1.0);
    EXPECT_NEAR(mid.x, 1.0, 1e-12);
    EXPECT_NEAR(mid.y, 0.5, 1e-12);
    EXPECT_NEAR(mid.psi, 0.2, 1e-12);
    EXPECT_NEAR(mid.kappa, 0.1, 1e-12);
}

TEST(PathSample, InterpolatesHeadingAcrossWrap) {
    bpt::Path path;
    path.points.push_back(bpt::PathPoint{0.0, 0.0, 0.0, bpt::kPi - 0.1, 0.0});
    path.points.push_back(bpt::PathPoint{1.0, -1.0, 0.0, -bpt::kPi + 0.1, 0.0});
    const bpt::PathPoint p = path.sample(0.5);
    EXPECT_NEAR(std::abs(p.psi), bpt::kPi, 1e-12);
}

TEST(PathSample, ClampsOutsideRange) {
    const bpt::Path path = two_point_path();
    const bpt::PathPoint before = path.sample(-3.0);
    const bpt::PathPoint after = path.sample(10.0);
    EXPECT_NEAR(before.x, 0.0, 1e-12);
    EXPECT_NEAR(before.y, 0.0, 1e-12);
    EXPECT_NEAR(after.x, 2.0, 1e-12);
    EXPECT_NEAR(after.y, 1.0, 1e-12);
}

TEST(PathSample, StraightLineMatchesArcLength) {
    const bpt::Path path = bpt::make_straight_line(10.0, 1.0);
    const bpt::PathPoint p = path.sample(2.5);
    EXPECT_NEAR(p.x, 2.5, 1e-12);
    EXPECT_NEAR(p.y, 0.0, 1e-12);
    EXPECT_NEAR(p.psi, 0.0, 1e-12);
}

// Path::project documents: "Signed cross-track error is positive to the left of the path tangent."
TEST(PathProject, CrossTrackSignOnHeadingNorthwest) {
    bpt::Path path;
    const double psi = 3.0 * bpt::kPi / 4.0;
    for (int i = 0; i <= 200; ++i) {
        const double s = 0.1 * static_cast<double>(i);
        path.points.push_back(bpt::PathPoint{s, s * std::cos(psi), s * std::sin(psi), psi, 0.0});
    }
    const double s = 10.0;
    const double bx = s * std::cos(psi);
    const double by = s * std::sin(psi);
    // Left normal of the tangent (cos psi, sin psi) is (-sin psi, cos psi).
    const double nx = -std::sin(psi);
    const double ny = std::cos(psi);

    const auto left = path.project(bx + 1.5 * nx, by + 1.5 * ny, -1.0, 10.0);
    const auto right = path.project(bx - 0.7 * nx, by - 0.7 * ny, -1.0, 10.0);
    EXPECT_NEAR(left.cte, 1.5, 1e-9);
    EXPECT_NEAR(left.s, s, 1e-9);
    EXPECT_NEAR(left.distance, 1.5, 1e-9);
    EXPECT_NEAR(right.cte, -0.7, 1e-9);
    EXPECT_NEAR(right.s, s, 1e-9);
}

TEST(PathProject, CrossTrackSignOnLeftTurningCircle) {
    const double radius = 20.0;
    const bpt::Path path = bpt::make_circle(radius, bpt::kPi / 2.0);
    // make_circle turns left around the centre (0, radius); inside the turn is the left side.
    const double theta = bpt::kPi / 4.0;
    const double cx = 0.0;
    const double cy = radius;
    const double ux = std::sin(theta);
    const double uy = -std::cos(theta);

    const auto inside = path.project(cx + (radius - 1.0) * ux, cy + (radius - 1.0) * uy, -1.0, 10.0);
    const auto outside = path.project(cx + (radius + 1.5) * ux, cy + (radius + 1.5) * uy, -1.0, 10.0);
    EXPECT_GT(inside.cte, 0.0);
    EXPECT_NEAR(inside.cte, 1.0, 1e-3);
    EXPECT_LT(outside.cte, 0.0);
    EXPECT_NEAR(outside.cte, -1.5, 1e-3);
    EXPECT_NEAR(inside.s, radius * theta, 1e-2);
    EXPECT_NEAR(inside.psi, theta, 1e-3);
}

TEST(PathProject, LocalWindowMatchesGlobalSearch) {
    const bpt::Path path = bpt::make_straight_to_curve();
    const bpt::PathPoint ref = path.sample(70.0);
    const double x = ref.x - 0.8 * std::sin(ref.psi);
    const double y = ref.y + 0.8 * std::cos(ref.psi);
    const auto global = path.project(x, y, -1.0, 10.0);
    const auto local = path.project(x, y, 68.0, 10.0);
    EXPECT_NEAR(global.s, 70.0, 1e-2);
    EXPECT_NEAR(local.s, global.s, 1e-9);
    EXPECT_NEAR(local.cte, global.cte, 1e-9);
    EXPECT_NEAR(global.cte, 0.8, 1e-3);
}

TEST(CurvePath, ArcCurvatureIsInverseRadius) {
    bpt::CurvePathSpec spec;
    spec.straight_lead = 30.0;
    spec.transition = 10.0;
    spec.radius = 15.0;
    spec.angle = 1.0;
    spec.straight_tail = 30.0;
    spec.sample_ds = 0.05;
    const bpt::Path path = bpt::make_curve_path(spec, "test_curve");

    const double arc_start = spec.straight_lead + spec.transition;
    const double arc_end = arc_start + spec.radius * spec.angle;
    const double margin = 0.5;
    int checked = 0;
    for (std::size_t i = 1; i + 1 < path.points.size(); ++i) {
        const bpt::PathPoint& p = path.points[i];
        if (p.s < arc_start + margin || p.s > arc_end - margin) {
            continue;
        }
        EXPECT_NEAR(p.kappa, 1.0 / spec.radius, 1e-9);
        // The geometry itself must turn at 1/radius, not only the stored kappa field.
        const bpt::PathPoint& prev = path.points[i - 1];
        const bpt::PathPoint& next = path.points[i + 1];
        const double dpsi = bpt::wrap_to_pi(next.psi - prev.psi);
        EXPECT_NEAR(dpsi / (next.s - prev.s), 1.0 / spec.radius, 1e-6);
        ++checked;
    }
    EXPECT_GT(checked, 100);

    const bpt::PathPoint mid = path.sample(0.5 * (arc_start + arc_end));
    EXPECT_NEAR(mid.kappa, 1.0 / spec.radius, 1e-9);
    const auto proj = path.project(mid.x, mid.y, -1.0, 10.0);
    EXPECT_NEAR(proj.kappa, 1.0 / spec.radius, 1e-9);
    EXPECT_NEAR(proj.cte, 0.0, 1e-9);
}

TEST(CurvePath, StraightSectionsHaveZeroCurvature) {
    bpt::CurvePathSpec spec;
    const bpt::Path path = bpt::make_curve_path(spec, "test_curve");
    const double arc_end = spec.straight_lead + spec.transition + spec.radius * spec.angle;
    for (const bpt::PathPoint& p : path.points) {
        if (p.s < spec.straight_lead - 0.1 || p.s > arc_end + 0.1) {
            EXPECT_NEAR(p.kappa, 0.0, 1e-12);
        }
    }
}
