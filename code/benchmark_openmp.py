"""Measure C++ direct, serial Barnes-Hut and OpenMP Barnes-Hut on shared inputs."""

import argparse
import csv
import os
from pathlib import Path
import random
import statistics
import tempfile

from compare_backends import max_relative_errors, random_particles, require_close, run, write_particles


THREADS = (1, 2, 4, 8)
FIELDS = ("n", "particles", "mode", "threads", "trial", "seconds", "ke_rel", "pe_rel", "angular_rel")


def run_case(binaries, work, n, steps, mode, threads):
    serial, parallel = binaries
    binary = parallel if mode == "openmp" else serial
    env = dict(os.environ, OMP_DYNAMIC="FALSE", OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads))
    method = "direct" if mode == "direct" else "barnes_hut"
    return run((binary, "cc"), work, n, steps, method, env=env)


def check(binaries):
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory) / "code"
        work.mkdir()
        write_particles(work / "particles.bin", random_particles(10, box=2.0, max_radius=0.12))
        _, reference = run_case(binaries, work, 10, 4, "serial", 1)
        for threads in (1, 4):
            _, values = run_case(binaries, work, 10, 4, "openmp", threads)
            require_close(values, reference)
    print("OpenMP diagnostics agree with serial Barnes-Hut.")


def measure(binaries, args):
    order = random.Random(2026)
    with tempfile.TemporaryDirectory() as directory, args.csv.open("w", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(FIELDS)
        work = Path(directory) / "code"
        work.mkdir()
        for n in args.exponents:
            write_particles(work / "particles.bin", random_particles(n))
            cases = [("serial", 1), *(("openmp", threads) for threads in THREADS)]
            if n <= args.direct_max:
                cases.append(("direct", 1))
            for trial in range(1, args.repeats + 1):
                order.shuffle(cases)
                samples = {case: run_case(binaries, work, n, args.steps, *case) for case in cases}
                reference = samples["serial", 1][1]
                for (mode, threads), (seconds, values) in samples.items():
                    if mode == "openmp":
                        require_close(values, reference)
                    errors = max_relative_errors(values, reference, args.steps)
                    writer.writerow((n, 1 << n, mode, threads, trial, seconds, *errors))
                output.flush()
            print(f"Measured N={1 << n:,}", flush=True)


def draw(input_csv, destination):
    import matplotlib.pyplot as plt

    plt.rcParams.update({"pdf.fonttype": 42, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    with input_csv.open(newline="") as source:
        groups = {}
        for row in csv.DictReader(source):
            key = int(row["particles"]), row["mode"], int(row["threads"])
            groups.setdefault(key, []).append(float(row["seconds"]))
    medians = {key: statistics.median(values) for key, values in groups.items()}
    sizes = sorted({n for n, _, _ in medians})
    destination.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.6, 3.5), layout="constrained")
    for mode, threads, label, color, marker in (
        ("direct", 1, "Direct pairwise", "#454950", "s"),
        ("serial", 1, "Barnes–Hut (serial)", "#17669A", "o"),
        ("openmp", 4, "Barnes–Hut + OpenMP (4 threads)", "#AB4B14", "^"),
    ):
        points = [(n, medians[n, mode, threads]) for n in sizes if (n, mode, threads) in medians]
        if points:
            ax.plot(*zip(*points), label=label, color=color, marker=marker, linewidth=1.8)
    ax.set(xscale="log", yscale="log", xlabel="Particle count", ylabel="Median wall time (s)")
    ax.set_xticks(sizes, [f"{n:,}" for n in sizes])
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.grid(axis="y", alpha=0.24)
    ax.legend(loc="upper left", frameon=False)
    fig.savefig(destination / "nbody_solver_runtime.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.6, 3.5), layout="constrained")
    colors = ("#675A76", "#286B9A", "#B34D19")
    for n, color, marker in zip(sizes[-3:], colors, ("o", "s", "^")):
        speedups = [medians[n, "serial", 1] / medians[n, "openmp", t] for t in THREADS]
        ax.plot(THREADS, speedups, color=color, marker=marker, linewidth=1.8, label=f"N = {n:,}")
    ax.axhline(1, color="#646970", linewidth=1, linestyle="--")
    ax.set(xlabel="OpenMP threads", ylabel="Speedup over serial Barnes–Hut")
    ax.set_xticks(THREADS)
    ax.grid(axis="y", alpha=0.24)
    ax.legend(loc="upper left", frameon=False)
    fig.savefig(destination / "nbody_openmp_speedup.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="benchmark data file")
    parser.add_argument("--binaries", type=Path, nargs=2, metavar=("SERIAL", "OPENMP"))
    parser.add_argument("--plot-dir", type=Path, help="render PDF plots from the CSV")
    parser.add_argument("--exponents", type=int, nargs="+", default=[10, 12, 14, 16, 17])
    parser.add_argument("--direct-max", type=int, default=14)
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not (args.binaries or args.plot_dir) or (args.check and (not args.binaries or args.plot_dir)):
        parser.error("supply --binaries to benchmark or check, or --plot-dir to render plots")
    if not args.check and args.csv is None:
        parser.error("supply --csv to benchmark or render plots")
    if args.steps < 1 or args.repeats < 1 or any(n < 0 or n > 30 for n in args.exponents):
        parser.error("steps and repeats must be positive and exponents must be between 0 and 30")
    if args.binaries:
        binaries = tuple(path.resolve() for path in args.binaries)
        if args.check:
            check(binaries)
        else:
            measure(binaries, args)
    if args.plot_dir:
        draw(args.csv, args.plot_dir)


if __name__ == "__main__":
    main()
