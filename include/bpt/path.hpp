#pragma once

#include <string>
#include <vector>

namespace bpt {

struct PathPoint {
    double s{0.0};
    double x{0.0};
    double y{0.0};
    double psi{0.0};
    double kappa{0.0};
};

struct CurvePathSpec {
    double straight_lead{40.0};
    double transition{12.0};
    double radius{20.0};
    double angle{1.15};
    double straight_tail{40.0};
    double sample_ds{0.05};
};

class Path {
public:
    std::string name;
    // Arc length of the designed geometry. Samples may continue past this
    // so a preview horizon still has a reference near the end.
    double nominal_length{0.0};
    std::vector<PathPoint> points;

    double length() const { return points.empty() ? 0.0 : points.back().s; }

    PathPoint sample(double s) const;

    struct Projection {
        double s{0.0};
        double x{0.0};
        double y{0.0};
        double psi{0.0};
        double kappa{0.0};
        double cte{0.0};
        double distance{0.0};
    };

    // Signed cross-track error is positive to the left of the path tangent.
    // A negative s_hint requests a search over the whole polyline.
    Projection project(double x, double y, double s_hint, double window) const;
};

Path make_straight_line(double length, double ds = 0.05);
Path make_circle(double radius, double angle, double ds = 0.05);
Path make_curve_path(const CurvePathSpec& spec, const std::string& name);
Path make_straight_to_curve();
Path make_figure_eight();
Path make_double_lane_change();

}  // namespace bpt
