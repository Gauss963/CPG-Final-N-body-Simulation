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
   - Try to optimize calculations using techniques like the Barnes-Hut algorithm for larger N to reduce computational complexity from O(N²) to O(N log N).

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
| 8    | 10/21  | Choose a topic before 10/28                          |
| 9    | 10/28  | Gathering papers on `N-body Simulation`              |
| 10   | 11/4   | Start initial coding, construct necessary functions. |
| 11   | 11/11  | Get ready to run the simulation                      |
| 12   | 11/18  | Run the simulation                                   |
| 13   | 11/25  | Run the simulation                                   |
| 14   | 12/2   | Finish plotting and statistics                       |
| 15   | 12/9   | Flash talk recording                                 |

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


## References

1. Hockney, R. W., & Eastwood, J. W. (1988). *Computer Simulation Using Particles*. CRC Press.
2. Aarseth, S. J. (2003). *Gravitational N-Body Simulations: Tools and Algorithms*. Cambridge University Press.
3. Binney, J., & Tremaine, S. (2008). *Galactic Dynamics* (2nd ed.). Princeton University Press.
4. Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.