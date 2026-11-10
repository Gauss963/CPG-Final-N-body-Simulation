# Computer Programming on Geosciences Final Project: `N-body Simulation`


## Topic
`N-body Simulation`

In physics and astronomy, an N-body simulation is a simulation of a dynamical system of particles, usually under the influence of physical forces, such as gravity. N-body simulations are widely used tools in astrophysics, from investigating the dynamics of few-body systems like the Earth-Moon-Sun system to understanding the evolution of the large-scale structure of the universe.


## Introduction and Expected Achievement
In this project, I aim to develop a computational model that simulates the gravitational interactions of a system of particles over time. The primary goal is to implement an N-body simulation that can accurately depict the dynamics of a system influenced by mutual gravitational forces. By the end of the project, I expect to:

- Develop a robust N-body simulation code using efficient numerical methods.
- Visualize the motion and interaction of particles within the system.
- Analyze the effects of varying initial conditions and parameters on system evolution.
- Gain deeper insights into complex gravitational systems and their behaviors.


## Name and Student ID
張德權 (GAUSS TE-CHUAN CHANG) B09501028


## Method
To achieve the project's goals, I will follow these steps:

1. **Programming Language and Tools**: Use `FORTRAN`, `C++` for its computations speed and `Matplotlib` for visualizations.

2. **Initialization**:
   - Set up initial conditions, including positions and velocities of particles.
   - Allow for different configurations, such as random distributions or specific arrangements like a disk galaxy.

3. **Force Calculation**:
   - Implement Newton's law of universal gravitation to calculate the force between each pair of particles.
   - The C++ and Fortran backends support Barnes–Hut gravity for larger N. They use a 3D octree to group distant masses, typically reducing the work to O(N log N).

4. **Time Integration**:
   - Use numerical integration methods such as the Verlet integration or Runge-Kutta methods to update particle positions and velocities over discrete time steps.

5. **Visualization**:
   - Create real-time animations to visualize the movement and interaction of particles.
   - Plot trajectories and possibly create 3D visualizations for more complex systems.

6. **Analysis**:
   - Examine energy conservation, angular momentum, and other physical quantities to verify the accuracy of the simulation.
   - Study phenomena such as orbit stability, collision rates, and clustering effects.


## Video link (youtube)

Simulation Animation on YouTube

https://www.youtube.com/watch?v=QOyjhbGaQ2A

[Simulation]
n = 7
L = 10.0
dt = 0.001
T = 10.0
G = 1.0

<!-- You are required to record a 3-min preview video to be available on 12/9 -->

## Planning and Timeline

| Week | Date   | Task                                                 |
| ---- | ------ | ---------------------------------------------------- |
| 9    | 10/28  | Gathering papers on `N-body Simulation`              |
| 10   | 11/4   | Start initial coding, construct necessary functions. |
| 11   | 11/11  | Get ready to run the simulation                      |
| 12   | 11/18  | Run the simulation                                   |
| 13   | 11/25  | Run the simulation                                   |
| 14   | 12/02  | Run the simulation                                   |
| 15   | 12/09  | Finish plotting and statistics                       |
| 16   | 12/16  | Flash all works                                      |

## User Guide

### Frontend

There is only 1 frontend in this project, the `Python Frontend` is used to do all the plots. The reason I did not use `PGPLOT` was because using `Matplotlib` was way easier and I can focus optimizing the backends.

### Backend

There are 3 backend you can use in this project. The `Python Backend`, `FORTRAN Backend`, and `C++ Backend`. The backends did all the simulation works and save all the datas in the `./data/*.bin`.

Simple description of how to run your simulation and analysis programs.

1. Compile `./code/Simulation.cc`, `./code/Simulation.f90`
2. run `./code/Simulation_cc && python3 ./code/PlotFunctions.py` or `./code/Simulation_FF && python3 ./code/PlotFunctions.py`
3. Check statistics plots in `./code/plot/`
4. run `python3 ./code/Simulation.py` and set `make_animation = True` and `python3 ./code/FolderActions.py` to make animation.
5. Check animation in `./animation/*.mp4`
6. Modify the parameter in `./code/config.ini`

#### Barnes–Hut in C++ and Fortran

For a short trial with 4,096 particles, set `n = 12`, `dt = 0.001` and `T = 0.01` in `code/config.ini`. The supplied timing values request 100 million steps. Build from the repository root and run the executable of your choice from `code/` so it can read that file.

```sh
clang++ -std=c++17 -O2 code/Simulation.cc code/random_mt.cc -o /tmp/cpg-nbody-cpp
gfortran -O2 code/Simulation.f90 -o /tmp/cpg-nbody-fortran
(cd code && /tmp/cpg-nbody-cpp)
(cd code && /tmp/cpg-nbody-fortran)
```

The supplied `code/config.ini` selects `gravity = barnes_hut` with `theta = 0.5` for both native backends. Smaller `theta` values open more cells, improving accuracy at the cost of speed. Set `theta = 0` to traverse the tree without approximating gravity, or `gravity = direct` to use pairwise gravity and collisions. Older configs without `gravity` use `direct`.

In [Barnes–Hut](https://www.nature.com/articles/324446a0) mode, the octree is rebuilt as the particles move. Distant cells that satisfy `cell width / distance to center of mass < theta` approximate both gravity and potential energy; nearby particles are evaluated individually. Collision detection uses the same tree to find possible contacts and checks those pairs exactly. Dense systems with many contacts can still require quadratic work.

The backends write `data/cc_energy_data_N_<N>.bin` and `data/FF_energy_data_N_<N>.bin` in the existing format, with one block of doubles each for kinetic energy, potential energy, and angular momentum magnitude. To plot the data, select the corresponding `simulator` in `code/PlotFunctions.py` and run it from `code/`.

Particle counts are `N = 2^n`. The native executables accept `0 <= n <= 30` with the standard 32-bit integer builds. The practical limit depends on available memory, the number of time steps, and the density of contacts. The three diagnostic arrays use `24 * floor(T / dt)` bytes before accounting for particles and the tree.

Run the comparison script to time both backends and compare their diagnostics over the same steps.

```sh
python3 code/compare_backends.py /tmp/cpg-nbody-cpp /tmp/cpg-nbody-fortran
python3 code/compare_backends.py /tmp/cpg-nbody-cpp /tmp/cpg-nbody-fortran --tree-only --exponents 17 20 --steps 1 --repeats 1
```

The script supplies both programs with the same optional `initial_conditions` file in the `[Simulation]` section. The file contains `N` records of eight native-endian 64-bit floats in order `x, y, z, vx, vy, vz, mass, radius`. Relative paths are resolved from the simulation's working directory. The benchmark uses fixed initial conditions, `theta = 0.5`, and `dt = 2^-17`. Its timings include process startup and I/O.

On an Apple M2 with 8 GiB of RAM, Clang and GNU Fortran builds at `-O2` both completed a one-step run with 1,048,576 particles (`n = 20`). It took about 42 seconds in C++ and 55 seconds in Fortran. This is the largest case tested. At 16,384 particles over four steps, Barnes–Hut took 0.32 seconds in C++ and 0.52 seconds in Fortran, versus 3.15 and 6.14 seconds with direct gravity. Those four-step timings are medians of three runs.

#### OpenMP in C++

On macOS with Apple Clang, install Homebrew's `libomp` and build the OpenMP executable from the repository root.

```sh
brew install libomp
omp_prefix="$(brew --prefix libomp)"
clang++ -std=c++17 -O2 -Xpreprocessor -fopenmp -I"$omp_prefix/include" \
  code/Simulation.cc code/random_mt.cc -L"$omp_prefix/lib" -Wl,-rpath,"$omp_prefix/lib" \
  -lomp -o /tmp/cpg-nbody-openmp
(cd code && OMP_NUM_THREADS=4 /tmp/cpg-nbody-openmp)
```

The OpenMP build distributes the Barnes–Hut gravity and potential calculations across threads for $N \ge 1024$. The tree build, collisions, and time integration run serially. Use `gravity = barnes_hut` in `config.ini` and set `OMP_NUM_THREADS` to choose the team size.

The benchmark compares direct pairwise, serial Barnes–Hut, and Barnes–Hut with 1, 2, 4, and 8 OpenMP threads on the same initial conditions. Run it with the serial executable built above. Figure generation uses Matplotlib and writes vector PDFs.

```sh
python3 code/benchmark_openmp.py --binaries /tmp/cpg-nbody-cpp /tmp/cpg-nbody-openmp \
  --csv /tmp/nbody_openmp_benchmark.csv --plot-dir /tmp/nbody-plots
```

The measurements used in the report are in `plot/nbody_openmp_benchmark.csv`. See `Tex/main.pdf` for the runtime and thread scaling figures.

### Continuous Integration and Release

GitHub Actions builds both native backends on a macOS runner: the C++ backend with Clang and the Fortran backend with Flang. The workflow runs for pull requests targeting `main` and pushes to `main`. After the first successful build on `main`, it creates the `v1.0.0` GitHub Release; subsequent builds leave that release unchanged.


## References

1. Hockney, R. W., & Eastwood, J. W. (1988). *Computer Simulation Using Particles*. CRC Press.
2. Aarseth, S. J. (2003). *Gravitational N-Body Simulations: Tools and Algorithms*. Cambridge University Press.
3. Binney, J., & Tremaine, S. (2008). *Galactic Dynamics* (2nd ed.). Princeton University Press.
4. Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.
