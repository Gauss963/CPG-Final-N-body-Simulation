#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <array>
#include <string>
#include <sstream>
#include <map>
#include <stdexcept>

struct Config {
    int n;
    double L;
    double dt;
    double T;
    double G;
};

Config read_config(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open config file: " + filename);
    }

    std::string line;
    bool in_simulation_section = false;
    std::map<std::string, std::string> kv;

    while (std::getline(file, line)) {
        while (!line.empty() && (line.back() == ' ' || line.back() == '\t')) line.pop_back();
        while (!line.empty() && (line.front() == ' ' || line.front() == '\t')) line.erase(line.begin());

        if (line.empty() || line[0] == ';' || line[0] == '#') {
            continue;
        }

        if (line[0] == '[' && line.back() == ']') {
            std::string section = line.substr(1, line.size()-2);
            if (section == "Simulation") {
                in_simulation_section = true;
            } else {
                in_simulation_section = false;
            }
            continue;
        }

        if (in_simulation_section) {
            // key = value
            auto pos = line.find('=');
            if (pos == std::string::npos) continue;
            std::string key = line.substr(0, pos);
            std::string val = line.substr(pos+1);

            while(!key.empty() && (key.back()==' '||key.back()=='\t')) key.pop_back();
            while(!key.empty() && (key.front()==' '||key.front()=='\t')) key.erase(key.begin());

            while(!val.empty() && (val.back()==' '||val.back()=='\t')) val.pop_back();
            while(!val.empty() && (val.front()==' '||val.front()=='\t')) val.erase(val.begin());

            kv[key] = val;
        }
    }

    Config c;
    {
        // n (int)
        if (kv.find("n") == kv.end()) throw std::runtime_error("Missing n in config");
        c.n = std::stoi(kv["n"]);

        // L (double)
        if (kv.find("L") == kv.end()) throw std::runtime_error("Missing L in config");
        c.L = std::stod(kv["L"]);

        // dt (double)
        if (kv.find("dt") == kv.end()) throw std::runtime_error("Missing dt in config");
        c.dt = std::stod(kv["dt"]);

        // T (double)
        if (kv.find("T") == kv.end()) throw std::runtime_error("Missing T in config");
        c.T = std::stod(kv["T"]);

        // G (double)
        if (kv.find("G") == kv.end()) throw std::runtime_error("Missing G in config");
        c.G = std::stod(kv["G"]);
    }

    return c;
}

struct PointMass {
    double mass;
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    std::array<double, 3> acceleration;
    double radius;
    
    PointMass(double m,
              const std::array<double,3>& pos,
              const std::array<double,3>& vel,
              double r)
    : mass(m), position(pos), velocity(vel), radius(r)
    {
        acceleration = {0.0, 0.0, 0.0};
    }

    void apply_force(const std::array<double,3>& force) {
        for (int i=0; i<3; i++) {
            acceleration[i] = force[i] / mass;
        }
    }

    void update(double dt) {
        for (int i=0; i<3; i++) {
            velocity[i] += acceleration[i] * dt;
            position[i] += velocity[i] * dt;
        }
    }
};

static inline double norm(const std::array<double,3>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline double dot(const std::array<double,3>& a, const std::array<double,3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline std::array<double,3> cross(const std::array<double,3>& a, const std::array<double,3>& b) {
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
}

static inline std::array<double,3> operator-(const std::array<double,3>& a, const std::array<double,3>& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}

static inline std::array<double,3> operator+(const std::array<double,3>& a, const std::array<double,3>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

static inline std::array<double,3> operator*(double s, const std::array<double,3>& v) {
    return {s*v[0], s*v[1], s*v[2]};
}

static inline std::array<double,3> operator/(const std::array<double,3>& v, double s) {
    return {v[0]/s, v[1]/s, v[2]/s};
}

bool check_collision(const PointMass& a, const PointMass& b) {
    double dist = norm(a.position - b.position);
    return dist <= (a.radius + b.radius);
}

void handle_collision(PointMass& a, PointMass& b) {
    // Compute normal
    auto normal = a.position - b.position;
    double dist = norm(normal);
    if (dist == 0.0) {
        return;
    }
    normal = (1.0/dist)*normal;
    // Relative velocity
    auto relative_velocity = a.velocity - b.velocity;
    double velocity_along_normal = dot(relative_velocity, normal);

    if (velocity_along_normal > 0) {
        // Objects are separating
        return;
    }

    double m1 = a.mass;
    double m2 = b.mass;
    double impulse_magnitude = (2.0 * velocity_along_normal) / (1.0/m1 + 1.0/m2);
    auto impulse = impulse_magnitude * normal;

    for (int i=0; i<3; i++) {
        a.velocity[i] -= (impulse[i] / m1);
        b.velocity[i] += (impulse[i] / m2);
    }
}

int main() {
    // Parameters
    Config config = read_config("config.ini");
    double L = config.L;
    double dt = config.dt;
    double T = config.T;
    double G = config.G;

    int N = (1 << config.n);
    int steps = static_cast<int>(T/dt);


    // Random initialization
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> mass_dist(1.0, 2.0);
    std::uniform_real_distribution<double> pos_dist(-L/2.0, L/2.0);
    std::uniform_real_distribution<double> vel_dist(-1.0, 1.0);
    std::uniform_real_distribution<double> rad_dist(0.001, 0.05);

    std::vector<PointMass> particles;
    particles.reserve(N);
    for (int i=0; i<N; i++){
        std::array<double,3> pos = {pos_dist(gen), pos_dist(gen), pos_dist(gen)};
        std::array<double,3> vel = {vel_dist(gen), vel_dist(gen), vel_dist(gen)};
        double m = mass_dist(gen);
        double r = rad_dist(gen);
        particles.emplace_back(m, pos, vel, r);
    }

    std::vector<double> total_kinetic_energy(steps,0.0);
    std::vector<double> total_potential_energy(steps,0.0);
    std::vector<double> total_angular_momentum(steps,0.0);

    // Since it may be necessary to analyze the particle position of each time step later, it can be stored as needed.
    // Formate: [x1, y1, z1, vx1, vy1, vz1, m1, r1, x2, y2, z2, vx2...]
    std::stringstream pc_filename;
    pc_filename << "../cc_data/particle_data_N_" << N << ".bin";
    std::ofstream particle_data(pc_filename.str(), std::ios::binary);
    // For each particle, 3pos + 3vel + 1mass + 1radius = 8 double per particle
    // steps * N * 8 * sizeof(double)

    for (int step = 0; step < steps; step++) {
        double progress = (double)step / steps * 100;
        std::cout << "Step " << step << "/" << steps << " Progress: " << progress << "\n" ;

        // Calculate Pairwise Gravitational Force
        std::vector<std::array<double,3>> forces(N, {0.0,0.0,0.0});
        for (int i=0; i<N; i++) {
            for (int j=i+1; j<N; j++) {
                auto r_vec = particles[j].position - particles[i].position;
                double r = norm(r_vec);
                if (r==0.0) continue;
                double F_mag = G * particles[i].mass * particles[j].mass / (r*r);
                auto r_hat = (1.0/r)*r_vec;
                auto F_i = F_mag * r_hat;
                auto F_j = (-1.0)*F_i;
                for (int k=0; k<3; k++){
                    forces[i][k] += F_i[k];
                    forces[j][k] += F_j[k];
                }
            }
        }

        // Apply force
        for (int i=0; i<N; i++) {
            particles[i].apply_force(forces[i]);
        }

        // Check and handle ollision
        for (int i=0; i<N; i++) {
            for (int j=i+1; j<N; j++) {
                if (check_collision(particles[i], particles[j])) {
                    handle_collision(particles[i], particles[j]);
                }
            }
        }

        // update
        for (int i=0; i<N; i++) {
            particles[i].update(dt);
        }

        // Calculate Kinetic Energy
        double KE = 0.0;
        for (auto &p : particles) {
            double v2 = dot(p.velocity, p.velocity);
            KE += 0.5 * p.mass * v2;
        }

        // Calculate Potential Energy
        double PE = 0.0;
        for (int i=0; i<N; i++) {
            for (int j=i+1; j<N; j++) {
                auto r_vec = particles[j].position - particles[i].position;
                double r = norm(r_vec);
                if (r != 0.0) {
                    PE += -G * particles[i].mass * particles[j].mass / r;
                }
            }
        }

        // Calculate Angular Momentum
        std::array<double,3> L_total = {0.0,0.0,0.0};
        for (auto &p : particles) {
            // L = r x (m v)
            auto mv = p.mass * p.velocity;
            auto Lp = cross(p.position, mv);
            for (int k=0; k<3; k++)
                L_total[k] += Lp[k];
        }
        double L_mag = norm(L_total);

        total_kinetic_energy[step] = KE;
        total_potential_energy[step] = PE;
        total_angular_momentum[step] = L_mag;

        // Save particle data in every frame
        // Output: pos(x,y,z), vel(vx,vy,vz), mass, radius
        for (auto &p : particles) {
            double buffer[8] = {p.position[0], p.position[1], p.position[2],
                                p.velocity[0], p.velocity[1], p.velocity[2],
                                p.mass, p.radius};
            // particle_data.write(reinterpret_cast<char*>(buffer), sizeof(buffer));
        }

    }

    particle_data.close();

    {   
        std::stringstream eg_filename;
        eg_filename << "../data/cc_energy_data_N_" << N << ".bin";
        std::ofstream energy_data(eg_filename.str(), std::ios::binary);
        energy_data.write((char*)total_kinetic_energy.data(), steps*sizeof(double));
        energy_data.write((char*)total_potential_energy.data(), steps*sizeof(double));
        energy_data.write((char*)total_angular_momentum.data(), steps*sizeof(double));
        energy_data.close();
    }

    std::cout << "Simulation finished.\n";
    return 0;
}
