#include "bpt/path.hpp"

#include "bpt/types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace bpt {
namespace {

PathPoint lerp_point(const PathPoint& a, const PathPoint& b, double t) {
    PathPoint p;
    p.s = (1.0 - t) * a.s + t * b.s;
    p.x = (1.0 - t) * a.x + t * b.x;
    p.y = (1.0 - t) * a.y + t * b.y;
    p.psi = wrap_to_pi(a.psi + t * wrap_to_pi(b.psi - a.psi));
    p.kappa = (1.0 - t) * a.kappa + t * b.kappa;
    return p;
}

void extend_straight(Path& path, double distance, double ds) {
    if (path.points.empty() || distance <= 0.0) {
        return;
    }
    const PathPoint last = path.points.back();
    const int steps = std::max(1, static_cast<int>(std::ceil(distance / ds)));
    const double step = distance / static_cast<double>(steps);
    double s = last.s;
    double x = last.x;
    double y = last.y;
    for (int i = 0; i < steps; ++i) {
        s += step;
        x += step * std::cos(last.psi);
        y += step * std::sin(last.psi);
        path.points.push_back(PathPoint{s, x, y, last.psi, 0.0});
    }
}

double smoothstep(double u) {
    // Zero value and zero first derivative at both ends, rises from 0 to 1.
    return u - std::sin(2.0 * kPi * u) / (2.0 * kPi);
}

}  // namespace

PathPoint Path::sample(double s) const {
    PathPoint empty;
    if (points.empty()) {
        return empty;
    }
    if (s <= points.front().s) {
        return points.front();
    }
    if (s >= points.back().s) {
        return points.back();
    }
    const auto upper = std::lower_bound(
        points.begin(), points.end(), s, [](const PathPoint& p, double value) { return p.s < value; });
    const std::size_t idx = static_cast<std::size_t>(upper - points.begin());
    const PathPoint& b = points[idx];
    const PathPoint& a = points[idx - 1];
    const double span = b.s - a.s;
    const double t = (span > 1e-12) ? (s - a.s) / span : 0.0;
    return lerp_point(a, b, std::clamp(t, 0.0, 1.0));
}

Path::Projection Path::project(double x, double y, double s_hint, double window) const {
    Projection best;
    best.distance = 1e300;
    if (points.size() < 2) {
        return best;
    }

    auto search = [&](std::size_t i0, std::size_t i1, Projection& out) {
        i1 = std::min(i1, points.size() - 1);
        for (std::size_t i = i0; i < i1; ++i) {
            const PathPoint& a = points[i];
            const PathPoint& b = points[i + 1];
            const double vx = b.x - a.x;
            const double vy = b.y - a.y;
            const double len2 = vx * vx + vy * vy;
            double t = 0.0;
            if (len2 > 1e-12) {
                t = ((x - a.x) * vx + (y - a.y) * vy) / len2;
                t = std::clamp(t, 0.0, 1.0);
            }
            const double px = a.x + t * vx;
            const double py = a.y + t * vy;
            const double dx = x - px;
            const double dy = y - py;
            const double dist2 = dx * dx + dy * dy;
            if (dist2 < out.distance) {
                const PathPoint mid = lerp_point(a, b, t);
                out.distance = dist2;
                out.s = mid.s;
                out.x = mid.x;
                out.y = mid.y;
                out.psi = mid.psi;
                out.kappa = mid.kappa;
                out.cte = -std::sin(mid.psi) * dx + std::cos(mid.psi) * dy;
            }
        }
    };

    const bool global = s_hint < 0.0;
    if (!global) {
        const double lo = s_hint - window;
        const double hi = s_hint + window;
        const auto begin = std::lower_bound(
            points.begin(), points.end(), lo, [](const PathPoint& p, double value) { return p.s < value; });
        std::size_t i0 = static_cast<std::size_t>(begin - points.begin());
        if (i0 > 0) {
            --i0;
        }
        const auto end = std::upper_bound(
            points.begin(), points.end(), hi, [](double value, const PathPoint& p) { return value < p.s; });
        const std::size_t i1 = static_cast<std::size_t>(end - points.begin());
        search(i0, i1, best);
        best.distance = std::sqrt(best.distance);
        if (best.distance <= 4.0) {
            return best;
        }
    }

    Projection global_best;
    global_best.distance = 1e300;
    search(0, points.size() - 1, global_best);
    global_best.distance = std::sqrt(global_best.distance);
    return global_best;
}

Path make_straight_line(double length, double ds) {
    Path path;
    path.name = "straight_line";
    const int steps = std::max(2, static_cast<int>(std::ceil(length / ds)));
    path.points.reserve(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i) {
        const double s = length * static_cast<double>(i) / static_cast<double>(steps);
        path.points.push_back(PathPoint{s, s, 0.0, 0.0, 0.0});
    }
    path.nominal_length = path.length();
    return path;
}

Path make_circle(double radius, double angle, double ds) {
    Path path;
    path.name = "circle";
    const double length = std::abs(radius * angle);
    const int steps = std::max(2, static_cast<int>(std::ceil(length / ds)));
    const double kappa = 1.0 / radius;
    path.points.reserve(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i) {
        const double theta = angle * static_cast<double>(i) / static_cast<double>(steps);
        PathPoint p;
        p.s = std::abs(radius) * std::abs(theta);
        p.x = radius * std::sin(theta);
        p.y = radius * (1.0 - std::cos(theta));
        p.psi = theta;
        p.kappa = kappa;
        path.points.push_back(p);
    }
    path.nominal_length = path.length();
    return path;
}

Path make_curve_path(const CurvePathSpec& spec, const std::string& name) {
    Path path;
    path.name = name;
    const double kappa_max = 1.0 / spec.radius;
    const double arc = std::abs(spec.radius * spec.angle);
    const double total = spec.straight_lead + spec.transition + arc + spec.straight_tail;
    const int steps = std::max(2, static_cast<int>(std::ceil(total / spec.sample_ds)));
    const double ds = total / static_cast<double>(steps);

    auto kappa_of = [&](double s) {
        if (s < spec.straight_lead) {
            return 0.0;
        }
        if (s < spec.straight_lead + spec.transition) {
            const double u = (s - spec.straight_lead) / spec.transition;
            return kappa_max * smoothstep(u);
        }
        if (s < spec.straight_lead + spec.transition + arc) {
            return kappa_max;
        }
        return 0.0;
    };

    double x = 0.0;
    double y = 0.0;
    double psi = 0.0;
    path.points.reserve(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i) {
        const double s = ds * static_cast<double>(i);
        const double kappa = kappa_of(std::min(s, total));
        path.points.push_back(PathPoint{s, x, y, wrap_to_pi(psi), kappa});
        if (i == steps) {
            break;
        }
        const double psi_mid = psi + 0.5 * kappa * ds;
        x += std::cos(psi_mid) * ds;
        y += std::sin(psi_mid) * ds;
        psi += kappa * ds;
    }
    path.nominal_length = path.length();
    extend_straight(path, 60.0, spec.sample_ds);
    return path;
}

Path make_straight_to_curve() {
    CurvePathSpec spec;
    spec.straight_lead = 45.0;
    spec.transition = 14.0;
    spec.radius = 22.0;
    spec.angle = 1.20;
    spec.straight_tail = 45.0;
    spec.sample_ds = 0.05;
    return make_curve_path(spec, "straight_to_curve");
}

Path make_figure_eight() {
    Path path;
    path.name = "figure_eight";
    const double scale = 64.0;
    const int samples = 2400;
    std::vector<PathPoint> lap;
    lap.reserve(static_cast<std::size_t>(samples) + 1);
    double s = 0.0;
    double prev_x = 0.0;
    double prev_y = 0.0;
    for (int i = 0; i <= samples; ++i) {
        const double t = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(samples);
        const double x = scale * std::sin(t);
        const double y = scale * std::sin(t) * std::cos(t);
        const double xp = scale * std::cos(t);
        const double yp = scale * std::cos(2.0 * t);
        const double xpp = -scale * std::sin(t);
        const double ypp = -2.0 * scale * std::sin(2.0 * t);
        const double speed2 = xp * xp + yp * yp;
        const double speed = std::sqrt(speed2);
        const double kappa = (xp * ypp - yp * xpp) / (speed2 * speed);
        if (i > 0) {
            s += std::hypot(x - prev_x, y - prev_y);
        }
        PathPoint point;
        point.s = s;
        point.x = x;
        point.y = y;
        point.psi = std::atan2(yp, xp);
        point.kappa = kappa;
        lap.push_back(point);
        prev_x = x;
        prev_y = y;
    }
    path.nominal_length = lap.back().s;
    // Three laps keep the preview on the loop after one full traversal.
    path.points.reserve(lap.size() * 3);
    for (int lap_index = 0; lap_index < 3; ++lap_index) {
        for (std::size_t i = 0; i < lap.size(); ++i) {
            if (lap_index > 0 && i == 0) {
                continue;
            }
            PathPoint point = lap[i];
            point.s += path.nominal_length * static_cast<double>(lap_index);
            path.points.push_back(point);
        }
    }
    return path;
}

Path make_double_lane_change() {
    Path path;
    path.name = "double_lane_change";
    const double ds = 0.05;
    const double x_start = 35.0;
    const double ramp = 18.0;
    const double hold = 28.0;
    const double lane = 3.6;
    const double x_end = x_start + 2.0 * ramp + hold + 50.0;

    auto blend = [](double u) { return smoothstep(std::clamp(u, 0.0, 1.0)); };
    auto blend_d = [](double u) {
        if (u <= 0.0 || u >= 1.0) {
            return 0.0;
        }
        return 1.0 - std::cos(2.0 * kPi * u);
    };
    auto blend_dd = [](double u) {
        if (u <= 0.0 || u >= 1.0) {
            return 0.0;
        }
        return 2.0 * kPi * std::sin(2.0 * kPi * u);
    };

    auto profile = [&](double x, double& y, double& yp, double& ypp) {
        y = 0.0;
        yp = 0.0;
        ypp = 0.0;
        if (x < x_start) {
            return;
        }
        if (x < x_start + ramp) {
            const double u = (x - x_start) / ramp;
            y = lane * blend(u);
            yp = lane * blend_d(u) / ramp;
            ypp = lane * blend_dd(u) / (ramp * ramp);
            return;
        }
        if (x < x_start + ramp + hold) {
            y = lane;
            return;
        }
        if (x < x_start + 2.0 * ramp + hold) {
            const double u = (x - x_start - ramp - hold) / ramp;
            y = lane * (1.0 - blend(u));
            yp = -lane * blend_d(u) / ramp;
            ypp = -lane * blend_dd(u) / (ramp * ramp);
            return;
        }
    };

    double s = 0.0;
    double prev_x = 0.0;
    double prev_y = 0.0;
    const int steps = static_cast<int>(std::ceil(x_end / ds));
    path.points.reserve(static_cast<std::size_t>(steps) + 1);
    for (int i = 0; i <= steps; ++i) {
        const double x = x_end * static_cast<double>(i) / static_cast<double>(steps);
        double y = 0.0;
        double yp = 0.0;
        double ypp = 0.0;
        profile(x, y, yp, ypp);
        if (i > 0) {
            s += std::hypot(x - prev_x, y - prev_y);
        }
        const double speed = std::sqrt(1.0 + yp * yp);
        PathPoint point;
        point.s = s;
        point.x = x;
        point.y = y;
        point.psi = std::atan2(yp, 1.0);
        point.kappa = ypp / (speed * speed * speed);
        path.points.push_back(point);
        prev_x = x;
        prev_y = y;
    }
    path.nominal_length = path.length();
    extend_straight(path, 50.0, ds);
    return path;
}

}  // namespace bpt
