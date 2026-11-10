#include "BarnesHut.h"

#include <iostream>
#include <random>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <map>
#include <stdexcept>

enum class GravityMethod { direct, barnes_hut };

struct Config {
    int n;
    double L;
    double dt;
    double T;
    double G;
    GravityMethod gravity = GravityMethod::direct;
    double theta = 0.5;
    std::string initial_conditions;
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
        while (!line.empty() && (line.back() == ' ' || line.back() == '\t' || line.back() == '\r')) line.pop_back();
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

    if (const auto method = kv.find("gravity"); method != kv.end()) {
        if (method->second == "barnes_hut") {
            c.gravity = GravityMethod::barnes_hut;
        } else if (method->second != "direct") {
            throw std::runtime_error("gravity must be direct or barnes_hut");
        }
    }
    if (const auto theta = kv.find("theta"); theta != kv.end()) {
        std::size_t parsed = 0;
        c.theta = std::stod(theta->second, &parsed);
        if (parsed != theta->second.size() || !std::isfinite(c.theta) || c.theta < 0.0) {
            throw std::runtime_error("theta must be a finite, non-negative number");
        }
    }
    if (const auto initial = kv.find("initial_conditions"); initial != kv.end()) {
        c.initial_conditions = initial->second;
    }
    return c;
}

static void initialize_particles(std::vector<PointMass>& particles, int count, const Config& config) {
    if (!config.initial_conditions.empty()) {
        std::ifstream input(config.initial_conditions, std::ios::binary | std::ios::ate);
        const auto bytes = static_cast<std::streamoff>(count) * 8 * sizeof(double);
        if (!input || input.tellg() != static_cast<std::streamoff>(bytes)) {
            throw std::runtime_error("Expected N records of 8 doubles in " + config.initial_conditions);
        }
        input.seekg(0);
        for (int i = 0; i < count; ++i) {
            std::array<double, 8> record;
            if (!input.read(reinterpret_cast<char*>(record.data()), sizeof(record)) ||
                !std::all_of(record.begin(), record.end(), [](double v) { return std::isfinite(v); }) ||
                record[6] <= 0 || record[7] < 0) {
                throw std::runtime_error("Invalid particle in " + config.initial_conditions);
            }
            particles.emplace_back(record[6], Vec3{record[0], record[1], record[2]},
                                   Vec3{record[3], record[4], record[5]}, record[7]);
        }
        return;
    }

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> mass_dist(1.0, 2.0);
    std::uniform_real_distribution<double> pos_dist(-config.L / 2.0, config.L / 2.0);
    std::uniform_real_distribution<double> vel_dist(-1.0, 1.0);
    std::uniform_real_distribution<double> rad_dist(0.001, 0.05);
    for (int i = 0; i < count; ++i) {
        Vec3 pos = {pos_dist(gen), pos_dist(gen), pos_dist(gen)};
        Vec3 vel = {vel_dist(gen), vel_dist(gen), vel_dist(gen)};
        double m = mass_dist(gen);
        double r = rad_dist(gen);
        particles.emplace_back(m, pos, vel, r);
    }
}

static Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0]};
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

int main() try {
    // Parameters
    Config config = read_config("config.ini");
    double L = config.L;
    double dt = config.dt;
    double T = config.T;
    double G = config.G;

    if (config.n < 0 || config.n >= std::numeric_limits<int>::digits) {
        throw std::runtime_error("n is outside the supported range for N = 2^n");
    }
    if (!std::isfinite(L) || L <= 0.0 || !std::isfinite(G)) {
        throw std::runtime_error("L must be positive and G must be finite");
    }
    if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(T) || T < 0.0 ||
        T / dt > std::numeric_limits<int>::max()) {
        throw std::runtime_error("Invalid dt or T, or too many time steps");
    }
    const int N = 1 << config.n;
    const int steps = static_cast<int>(T / dt);
    const bool use_tree = config.gravity == GravityMethod::barnes_hut;

    std::vector<PointMass> particles;
    particles.reserve(N);
    initialize_particles(particles, N, config);

    std::vector<double> total_kinetic_energy(steps,0.0);
    std::vector<double> total_potential_energy(steps,0.0);
    std::vector<double> total_angular_momentum(steps,0.0);

    BarnesHut tree(particles);
    std::vector<std::size_t> collisions;
    const auto evaluate_gravity = [&] {
        if (!use_tree) return direct_gravity(particles, G);
        tree.rebuild();
        return tree.gravity(G, config.theta);
    };

    std::cout << "Gravity: " << (use_tree ? "barnes_hut" : "direct");
    if (use_tree) std::cout << " (theta=" << config.theta << ')';
    std::cout << '\n';

    // Accelerations at the initial positions; subsequent evaluations also record energy.
    if (steps > 0) evaluate_gravity();

    for (int step = 0; step < steps; step++) {
        double progress = (double)step / steps * 100;
        std::cout << "Step " << step << "/" << steps << " Progress: " << progress << "\n" ;

        // Resolve contacts at the current positions before updating velocities and positions.
        for (std::size_t i = 0; i < particles.size(); ++i) {
            if (use_tree) {
                tree.collisions_after(i, collisions);
                for (const auto j : collisions) handle_collision(particles[i], particles[j]);
            } else {
                for (std::size_t j = i + 1; j < particles.size(); ++j) {
                    if (check_collision(particles[i], particles[j])) {
                        handle_collision(particles[i], particles[j]);
                    }
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

        // At the new positions, compute potential energy and prepare the next step.
        const double PE = evaluate_gravity();

        // Calculate Angular Momentum
        Vec3 L_total{};
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
    }

    std::filesystem::create_directories("../data");
    const auto filename = "../data/cc_energy_data_N_" + std::to_string(N) + ".bin";
    std::ofstream energy_data(filename, std::ios::binary);
    if (!energy_data) throw std::runtime_error("Could not open output file: " + filename);
    const std::streamsize bytes = static_cast<std::streamsize>(steps) * sizeof(double);
    if (steps > 0) {
        energy_data.write(reinterpret_cast<const char*>(total_kinetic_energy.data()), bytes);
        energy_data.write(reinterpret_cast<const char*>(total_potential_energy.data()), bytes);
        energy_data.write(reinterpret_cast<const char*>(total_angular_momentum.data()), bytes);
    }
    energy_data.close();
    if (!energy_data) throw std::runtime_error("Could not write output file: " + filename);

    std::cout << "Simulation finished.\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "Simulation error: " << error.what() << '\n';
    return 1;
}
