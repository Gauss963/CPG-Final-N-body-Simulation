"""Run the native solvers on identical initial conditions and compare diagnostics."""

import argparse
import math
from pathlib import Path
import random
import statistics
import struct
import subprocess
import tempfile
import time


PARTICLE = struct.Struct("=8d")
DT = 2**-17


def random_particles(n, box=10.0, max_radius=0.05):
    rng = random.Random(42)
    for _ in range(1 << n):
        pos = [rng.uniform(-box / 2, box / 2) for _ in range(3)]
        vel = [rng.uniform(-1, 1) for _ in range(3)]
        yield (*pos, *vel, rng.uniform(1, 2), rng.uniform(0.001, max_radius))


def write_particles(path, particles):
    with path.open("wb") as output:
        for particle in particles:
            output.write(PARTICLE.pack(*particle))


def run(backend, work, n, steps, method, theta=0.5, dt=DT, env=None):
    binary, prefix = backend
    (work / "config.ini").write_text(
        f"[Simulation]\nn = {n}\nL = 10\ndt = {dt}\nT = {dt * steps}\nG = 1\n"
        f"gravity = {method}\ntheta = {theta}\ninitial_conditions = particles.bin\n"
    )
    start = time.perf_counter()
    subprocess.run([binary], cwd=work, env=env, stdout=subprocess.DEVNULL, check=True)
    seconds = time.perf_counter() - start
    output = work.parent / "data" / f"{prefix}_energy_data_N_{1 << n}.bin"
    values = struct.unpack(f"={3 * steps}d", output.read_bytes())
    if not all(map(math.isfinite, values)):
        raise ValueError(f"Non-finite diagnostics from {binary} ({method})")
    return seconds, values


def max_relative_errors(values, reference, steps):
    return tuple(
        max(abs(a - b) / max(abs(b), 1e-30) for a, b in zip(values[i:i + steps], reference[i:i + steps]))
        for i in range(0, 3 * steps, steps)
    )


def require_close(actual, expected):
    for a, b in zip(actual, expected):
        if not math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-12):
            raise ValueError(f"Diagnostics disagree: got {a}, expected {b}")


def check(backends, work):
    particles = work / "particles.bin"
    # Two touching bodies approach with speeds 1. After bouncing and applying
    # gravity, their separation is 2 + 2 * dt * (1 - dt / 4).
    write_particles(particles, [(-1, 0, 0, 1, 0, 0, 1, 1), (1, 0, 0, -1, 0, 0, 1, 1)])
    dt = 0.125
    speed = 1 - dt / 4
    expected = (speed**2, -1 / (2 + 2 * dt * speed), 0.0)
    for backend in backends:
        for method in ("direct", "barnes_hut"):
            require_close(run(backend, work, 1, 1, method, dt=dt)[1], expected)

    # This dense fixture exercises collision pruning and impulse order across cells.
    write_particles(particles, random_particles(8, box=1.0, max_radius=0.25))
    _, exact = run(backends[0], work, 8, 4, "direct")
    _, approx = run(backends[0], work, 8, 4, "barnes_hut")
    require_close(run(backends[1], work, 8, 4, "direct")[1], exact)
    require_close(run(backends[1], work, 8, 4, "barnes_hut")[1], approx)
    for backend in backends:
        require_close(run(backend, work, 8, 4, "barnes_hut", theta=0)[1], exact)
    print("Both backends passed the analytic collision and shared initial state checks.")


def benchmark(backends, work, args):
    methods = ("barnes_hut",) if args.tree_only else ("direct", "barnes_hut")
    for n in args.exponents:
        write_particles(work / "particles.bin", random_particles(n))
        timings, diagnostics = {}, {}
        for _ in range(args.repeats):
            for method in methods:
                for backend in backends:
                    seconds, values = run(backend, work, n, args.steps, method)
                    key = backend[1], method
                    timings.setdefault(key, []).append(seconds)
                    diagnostics[key] = values

        print(f"N={1 << n:,} ({args.steps} steps, median of {args.repeats} runs)")
        for method in methods:
            cpp = statistics.median(timings["cc", method])
            ff = statistics.median(timings["FF", method])
            print(f"  {method:10s} C++ {cpp:.4f} s, Fortran {ff:.4f} s")
            errors = max_relative_errors(diagnostics["FF", method], diagnostics["cc", method], args.steps)
            print(f"    Max Fortran vs C++ relative difference (KE, PE, |L|): {errors[0]:.2e}, {errors[1]:.2e}, {errors[2]:.2e}")
        if not args.tree_only:
            for name in ("cc", "FF"):
                errors = max_relative_errors(diagnostics[name, "barnes_hut"], diagnostics[name, "direct"], args.steps)
                print(f"  {name} max tree vs direct relative difference (KE, PE, |L|): "
                      f"{errors[0]:.2e}, {errors[1]:.2e}, {errors[2]:.2e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cpp", type=Path, help="compiled C++ executable")
    parser.add_argument("fortran", type=Path, help="compiled Fortran executable")
    parser.add_argument("--exponents", type=int, nargs="+", default=[10, 12, 14], help="particle counts as powers of two")
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--tree-only", action="store_true", help="skip the quadratic solver at large N")
    parser.add_argument("--check", action="store_true", help="run numerical regression checks instead of timing")
    args = parser.parse_args()
    if args.steps < 1 or args.repeats < 1 or any(n < 0 or n > 30 for n in args.exponents):
        parser.error("steps and repeats must be positive, and exponents must be between 0 and 30")
    backends = ((args.cpp.resolve(), "cc"), (args.fortran.resolve(), "FF"))
    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory) / "code"
        work.mkdir()
        if args.check:
            check(backends, work)
        else:
            benchmark(backends, work, args)


if __name__ == "__main__":
    main()
