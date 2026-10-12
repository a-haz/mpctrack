#!/usr/bin/env python3
"""Draw simulated path-tracking figures from the benchmark CSV files."""

import argparse
import csv
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MPC_COLOR = "#0072B2"
PP_COLOR = "#D55E00"
PATH_COLOR = "#222222"


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", default="docs/benchmark.csv")
    parser.add_argument("--trajectories", default="docs/trajectories.csv")
    parser.add_argument("--paths", default="docs/paths.csv")
    parser.add_argument("--outdir", default="docs")
    args = parser.parse_args()

    benchmark = read_csv(args.benchmark)
    trajectories = read_csv(args.trajectories)
    path_rows = read_csv(args.paths)

    paths = []
    for row in path_rows:
        if row["path"] not in paths:
            paths.append(row["path"])
    speeds = []
    for row in benchmark:
        if row["speed_mps"] not in speeds:
            speeds.append(row["speed_mps"])

    path_xy = defaultdict(lambda: {"x": [], "y": []})
    for row in path_rows:
        path_xy[row["path"]]["x"].append(float(row["x"]))
        path_xy[row["path"]]["y"].append(float(row["y"]))

    grouped = defaultdict(lambda: {"x": [], "y": [], "t": [], "cte": []})
    for row in trajectories:
        key = (row["controller"], row["path"], row["speed_mps"])
        grouped[key]["x"].append(float(row["x"]))
        grouped[key]["y"].append(float(row["y"]))
        grouped[key]["t"].append(float(row["t"]))
        grouped[key]["cte"].append(float(row["cte"]))

    overlay_speed = speeds[len(speeds) // 2]
    fig, axes = plt.subplots(len(paths), 1, figsize=(11.0, 11.5), constrained_layout=True)
    if len(paths) == 1:
        axes = [axes]
    for axis, path_name in zip(axes, paths):
        axis.plot(path_xy[path_name]["x"], path_xy[path_name]["y"], color=PATH_COLOR, linewidth=1.6,
                  linestyle="--", label="reference", zorder=1)
        for controller, color, label in (
            ("ltv_mpc", MPC_COLOR, "LTV-MPC"),
            ("pure_pursuit", PP_COLOR, "pure pursuit"),
        ):
            series = grouped[(controller, path_name, overlay_speed)]
            axis.plot(series["x"], series["y"], color=color, linewidth=1.4, label=label, zorder=2)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(f"{path_name.replace('_', ' ')} at {overlay_speed} m/s")
        axis.set_xlabel("x (m)")
        axis.set_ylabel("y (m)")
        axis.grid(True, alpha=0.3)
        axis.legend(loc="best", frameon=False)
    fig.suptitle("Simulated paths", fontsize=14)
    fig.savefig(f"{args.outdir}/paths_overlay.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(len(paths), len(speeds), figsize=(12.0, 8.2), sharey=True, constrained_layout=True)
    if len(paths) == 1 and len(speeds) == 1:
        axes = [[axes]]
    elif len(paths) == 1:
        axes = [axes]
    elif len(speeds) == 1:
        axes = [[axis] for axis in axes]
    for row, path_name in enumerate(paths):
        for col, speed in enumerate(speeds):
            axis = axes[row][col]
            for controller, color, label in (
                ("ltv_mpc", MPC_COLOR, "LTV-MPC"),
                ("pure_pursuit", PP_COLOR, "pure pursuit"),
            ):
                series = grouped[(controller, path_name, speed)]
                axis.plot(series["t"], series["cte"], color=color, linewidth=1.2, label=label)
            axis.axhline(0.0, color="#666666", linewidth=0.6)
            axis.grid(True, alpha=0.3)
            if row == 0:
                axis.set_title(f"{speed} m/s")
            if col == 0:
                axis.set_ylabel(f"{path_name.replace('_', ' ')}\ncross-track (m)")
            if row == len(paths) - 1:
                axis.set_xlabel("time (s)")
            if row == 0 and col == 0:
                axis.legend(loc="best", frameon=False)
    fig.suptitle("Simulated cross-track error", fontsize=14)
    fig.savefig(f"{args.outdir}/cte_over_time.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, len(paths), figsize=(12.0, 3.8), sharey=True, constrained_layout=True)
    if len(paths) == 1:
        axes = [axes]
    width = 0.36
    for axis, path_name in zip(axes, paths):
        mpc = []
        pursuit = []
        for speed in speeds:
            mpc_row = next(row for row in benchmark
                           if row["controller"] == "ltv_mpc" and row["path"] == path_name
                           and row["speed_mps"] == speed)
            pp_row = next(row for row in benchmark
                          if row["controller"] == "pure_pursuit" and row["path"] == path_name
                          and row["speed_mps"] == speed)
            mpc.append(float(mpc_row["rms_cte_m"]))
            pursuit.append(float(pp_row["rms_cte_m"]))
        positions = list(range(len(speeds)))
        axis.bar([p - width / 2 for p in positions], mpc, width=width, color=MPC_COLOR, label="LTV-MPC")
        axis.bar([p + width / 2 for p in positions], pursuit, width=width, color=PP_COLOR, label="pure pursuit")
        axis.set_xticks(positions, [f"{speed}" for speed in speeds])
        axis.set_xlabel("speed (m/s)")
        axis.set_title(path_name.replace("_", " "))
        axis.grid(True, axis="y", alpha=0.3)
    axes[0].set_ylabel("RMS cross-track error (m)")
    axes[0].legend(frameon=False)
    fig.suptitle("Simulated RMS cross-track error", fontsize=14)
    fig.savefig(f"{args.outdir}/rms_cte.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
