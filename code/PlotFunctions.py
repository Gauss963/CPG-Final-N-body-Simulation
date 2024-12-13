import configparser
import matplotlib.pyplot as plt
import numpy as np

def plot_from_simulator():

    config = configparser.ConfigParser()
    config.read('config.ini')

    n = config.getint('Simulation', 'n')
    L = config.getfloat('Simulation', 'L')
    dt = config.getfloat('Simulation', 'dt')
    T = config.getfloat('Simulation', 'T')
    G = config.getfloat('Simulation', 'G')

    N = (1 << n)
    steps = int(T/dt)
    
    # simulator = 'py'
    # simulator = 'cc'
    simulator = 'FF'

    # Read from the binary
    with open(f'../data/{simulator}_energy_data_N_{N}.bin', "rb") as f:
        data = np.fromfile(f, dtype=np.float64, count=3*steps)
        total_kinetic_energy = data[0:steps]
        total_potential_energy = data[steps:2*steps]
        total_angular_momentum = data[2*steps:3*steps]

    time_array = np.arange(steps) * dt

    # Plot Energy v.s. Time Step
    plt.figure(figsize = (10, 5))
    plt.plot(time_array, total_kinetic_energy, label='Total Kinetic Energy', linestyle='--')
    plt.plot(time_array, total_potential_energy, label='Total Potential Energy', linestyle='--')
    plt.plot(time_array, total_kinetic_energy + total_potential_energy, label='Total Energy', linestyle='-.')
    plt.xlabel("Time (s)")
    plt.ylabel("Energy")
    plt.title(f'Energy vs Time, $N = {N}$')
    plt.grid(True)
    plt.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig(f'../plot/{simulator}_energy_vs_time_{N}_points.pdf')

    # Plot Angular Momentum v.s. Time Step
    plt.figure(figsize=(10,5))
    plt.plot(time_array, total_angular_momentum, label='|Total Angular Momentum|')
    plt.xlabel("Time (s)")
    plt.ylabel("Angular Momentum Magnitude")
    plt.title(f'Total Angular Momentum vs Time, $N = {N}$')
    plt.grid(True)
    plt.legend(loc = 'upper right')
    plt.tight_layout()
    plt.savefig(f'../plot/{simulator}_angular_momentum_vs_time_{N}_points.pdf')
    
    
if __name__ == "__main__":
    plot_from_simulator()