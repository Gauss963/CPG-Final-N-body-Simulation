import matplotlib.pyplot as plt
import numpy as np

import FolderActions
import PointMass


sim_seed = 42
np.random.seed(sim_seed)


make_animation = False

# Define Parameters
n = 7
N = 2**n      # Number of "PointMass"
L = 10.0      # Simulation Space Size (initial)
dt = 0.001    # Time Step Len
T = 10.0      # Simulation Life Time
steps = int(T / dt)

# Gravitational Constant
G = 1.0

# Generate N "PointMass" Randomly
masses     = np.random.uniform(1.0, 2.0, size = N)          # Mass Range
positions  = np.random.uniform(-L/2, L/2, size = (N, 3))    # Position Range
velocities = np.random.uniform(-1.0, 1.0, size = (N, 3))    # Velocitiy Range
radii      = np.random.uniform(0.001, 0.05, size = N)       # Radius Range

particles = [PointMass.PointMass(mass=masses[i],
                                 position=positions[i],
                                 velocity=velocities[i],
                                 radius=radii[i])
             for i in range(N)]


total_kinetic_energy = np.zeros(steps)
total_potential_energy = np.zeros(steps)
total_angular_momentum = np.zeros(steps)

for step in range(steps):
    
    print(f"Step {step}/{steps}")
    
    # Calculate Gravitational Force
    forces = [np.zeros(3) for _ in range(N)]
    
    # Calculate Pairwise Gravitational Force
    for i in range(N):
        for j in range(i+1, N):
            r_vec = particles[j].get_position() - particles[i].get_position()
            r = np.linalg.norm(r_vec)
            if r == 0:
                continue

            F_mag = G * particles[i].get_mass() * particles[j].get_mass() / (r**2)

            r_hat = r_vec / r

            F_i = F_mag * r_hat

            F_j = -F_i
            forces[i] += F_i
            forces[j] += F_j


    for i in range(N):
        particles[i].apply_force(forces[i])

    # Check And Handle Collision
    for i in range(N):
        for j in range(i+1, N):
            if particles[i].check_collision(particles[j]):
                particles[i].handle_collision(particles[j])

    # Update Position and Velocitiy
    for i in range(N):
        particles[i].update(dt)
        
    # Calculate Kinetic Energy
    KE = 0.0
    for p in particles:
        v = p.get_velocity()
        m = p.get_mass()
        KE += 0.5 * m * np.dot(v, v)

    # Calculate Potential Energy
    PE = 0.0
    for i in range(N):
        for j in range(i + 1, N):
            r_vec = particles[j].get_position() - particles[i].get_position()
            r = np.linalg.norm(r_vec)
            if r != 0:
                PE += -G * particles[i].get_mass() * particles[j].get_mass() / r
    
    # Calculate Angular Momentum
    # L = sum over i of (r_i x (m_i v_i))
    L_total = np.zeros(3)
    for p in particles:
        r = p.get_position()
        v = p.get_velocity()
        m = p.get_mass()
        L_total += np.cross(r, m*v)
    L_mag = np.linalg.norm(L_total)

    total_kinetic_energy[step] = KE
    total_potential_energy[step] = PE
    total_angular_momentum[step] = L_mag
    
    if make_animation:
        
        x_positions = [p.get_position()[0] for p in particles]
        y_positions = [p.get_position()[1] for p in particles]

        plt.figure(figsize = (8, 8))
        plt.scatter(x_positions, y_positions, s = 10, c = 'b', alpha = 0.5)
        plt.xlim(-L * 2, L * 2)
        plt.ylim(-L * 2, L * 2)
        plt.title(f"Step {step}/{steps}")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(True)
        
        plt.tight_layout()
        plt.savefig(f'../animation/frame/{step:04d}.png', dpi = 300)
        plt.close()
        


    if step % (steps // 10) == 0:
        print(f"Step {step}/{steps}")
        
        x_positions = [p.get_position()[0] for p in particles]
        y_positions = [p.get_position()[1] for p in particles]

        plt.figure(figsize = (8, 8))
        plt.scatter(x_positions, y_positions, s = 10, c='b', alpha = 0.5)
        plt.xlim(-L * 2, L * 2)
        plt.ylim(-L * 2, L * 2)
        plt.title(f"Step {step}/{steps}")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(True)
        
        plt.tight_layout()
        plt.savefig(f'../plot/step_{step}.pdf', dpi = 300)
        plt.close()



# Plot Energy v.s. Time Step
time_array = np.arange(steps) * dt

plt.figure(figsize = (10, 5))
plt.plot(time_array, total_kinetic_energy, label = 'Total Kinetic Energy', linestyle='--')
plt.plot(time_array, total_potential_energy, label = 'Total Potential Energy', linestyle='--')
plt.plot(time_array, total_kinetic_energy + total_potential_energy, label = 'Total Energy', linestyle='-.')
plt.xlabel("Time (s)")
plt.ylabel("Energy")
plt.title("Energy vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'../plot/py_energy_vs_time_{N}_points.pdf')

# Plot Angular Momentum v.s. Time Step
plt.figure(figsize = (10, 5))
plt.plot(time_array, total_angular_momentum, label = '|Total Angular Momentum|')
plt.xlabel("Time (s)")
plt.ylabel("Angular Momentum Magnitude")
plt.title("Total Angular Momentum vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(f'../plot/py_angular_momentum_vs_time_{N}_points.pdf')



if make_animation:
    FolderActions.make_animation(N)
FolderActions.delete_pycache()

print("Simulation finished.")