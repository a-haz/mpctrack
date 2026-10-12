#include "bpt/control.hpp"
#include "bpt/simulator.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string argument(int argc, char** argv, const std::string& key, const std::string& fallback) {
    for (int i = 1; i + 1 < argc; ++i) {
        if (key == argv[i]) {
            return argv[i + 1];
        }
    }
    return fallback;
}

unsigned scenario_seed(const std::string& name, double speed) {
    unsigned hash = 2166136261u;
    for (unsigned char character : name) {
        hash ^= character;
        hash *= 16777619u;
    }
    hash ^= static_cast<unsigned>(std::lround(speed * 10.0));
    hash *= 16777619u;
    return hash;
}

void write_paths(const std::string& filename, const std::vector<bpt::Path>& paths) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("could not write " + filename);
    }
    out << "path,s,x,y\n";
    out << std::fixed << std::setprecision(5);
    for (const bpt::Path& path : paths) {
        const int stride = 4;
        for (std::size_t i = 0; i < path.points.size(); i += static_cast<std::size_t>(stride)) {
            const bpt::PathPoint& point = path.points[i];
            if (point.s > path.nominal_length + 5.0) {
                break;
            }
            out << path.name << "," << point.s << "," << point.x << "," << point.y << "\n";
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    const std::string csv_path = argument(argc, argv, "--csv", "docs/benchmark.csv");
    const std::string trajectory_path = argument(argc, argv, "--trajectories", "docs/trajectories.csv");
    const std::string paths_path = argument(argc, argv, "--paths", "docs/paths.csv");

    const bpt::VehicleParams vehicle;
    const bpt::MpcParams mpc_tuning;
    const bpt::PurePursuitParams pursuit_tuning;
    const std::vector<bpt::Path> paths = {bpt::make_straight_to_curve(), bpt::make_figure_eight(),
                                          bpt::make_double_lane_change()};
    const std::vector<double> speeds = {6.0, 10.0, 14.0};
    const double s0 = 8.0;

    write_paths(paths_path, paths);

    std::ofstream csv(csv_path);
    std::ofstream trajectories(trajectory_path);
    if (!csv || !trajectories) {
        std::cerr << "Could not open benchmark outputs\n";
        return 1;
    }
    csv << "controller,path,speed_mps,rms_cte_m,max_cte_m,rms_heading_deg,max_heading_deg,"
           "mean_solve_ms,max_solve_ms,constraint_violations\n";
    trajectories << "controller,path,speed_mps,t,x,y,cte,heading_deg\n";
    trajectories << std::fixed << std::setprecision(6);

    std::cout << std::left << std::setw(14) << "controller" << std::setw(22) << "path" << std::setw(8) << "speed"
              << std::setw(12) << "rms_cte" << std::setw(12) << "max_cte" << std::setw(12) << "rms_head"
              << std::setw(12) << "solve_ms" << "violations\n";

    for (const bpt::Path& path : paths) {
        std::cout << path.name << " nominal length " << path.nominal_length << " m, stored " << path.length()
                  << " m\n";
        for (double speed : speeds) {
            const double duration = std::max(0.5, (path.nominal_length - s0) / speed);
            bpt::SimConfig config;
            config.dt = mpc_tuning.dt;
            config.duration = duration;
            config.s0 = s0;
            config.speed0 = speed;
            config.violation_tolerance = 1e-4;
            config.noise.enabled = true;
            config.noise.seed = scenario_seed(path.name, speed);

            for (int controller_index = 0; controller_index < 2; ++controller_index) {
                bpt::SimResult result;
                std::string name;
                if (controller_index == 0) {
                    bpt::LtvMpc controller(path, vehicle, mpc_tuning, speed);
                    name = controller.name();
                    result = bpt::simulate(controller, path, vehicle, config);
                } else {
                    bpt::PurePursuit controller(path, vehicle, pursuit_tuning, speed, config.dt);
                    name = controller.name();
                    result = bpt::simulate(controller, path, vehicle, config);
                }

                csv << name << "," << path.name << "," << std::fixed << std::setprecision(0) << speed << ","
                    << std::setprecision(6) << result.rms_cte << "," << result.max_cte << ","
                    << result.rms_heading_deg << "," << result.max_heading_deg << "," << std::setprecision(4)
                    << result.mean_solve_ms << "," << result.max_solve_ms << "," << result.violations << "\n";

                for (const bpt::Sample& sample : result.samples) {
                    trajectories << name << "," << path.name << "," << std::fixed << std::setprecision(0) << speed
                                 << "," << std::setprecision(6) << sample.t << "," << sample.x << "," << sample.y
                                 << "," << sample.cte << "," << sample.heading_deg << "\n";
                }

                std::cout << std::left << std::setw(14) << name << std::setw(22) << path.name << std::setw(8)
                          << speed << std::setw(12) << result.rms_cte << std::setw(12) << result.max_cte
                          << std::setw(12) << result.rms_heading_deg << std::setw(12) << result.mean_solve_ms
                          << result.violations << "\n";
            }
        }
    }
    return 0;
}
