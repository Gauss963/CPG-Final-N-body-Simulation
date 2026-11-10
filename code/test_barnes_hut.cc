#include "BarnesHut.h"

#include <iostream>
#include <random>
#include <stdexcept>

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static bool near(double a, double b) {
    return std::abs(a - b) <= 1e-11 * std::max({1.0, std::abs(a), std::abs(b)});
}

static void check_collisions(const BarnesHut& tree, const std::vector<PointMass>& particles) {
    std::vector<std::size_t> expected, actual;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        expected.clear();
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            if (check_collision(particles[i], particles[j])) expected.push_back(j);
        }
        tree.collisions_after(i, actual);
        require(actual == expected, "Tree missed, duplicated or reordered a collision");
    }
}

static void test_small_systems() {
    std::vector<PointMass> particles{{2, {0, 0, 0}, {}, 1}};
    BarnesHut tree(particles);
    tree.rebuild();
    require(tree.gravity(2, 0.5) == 0 && particles[0].acceleration == Vec3{},
            "One particle must have no self-gravity");

    particles.emplace_back(3, Vec3{2, 0, 0}, Vec3{}, 1);
    tree.rebuild();
    require(near(tree.gravity(2, 0.5), -6), "Incorrect two-body potential energy");
    require(near(particles[0].acceleration[0], 1.5) && near(particles[1].acceleration[0], -1),
            "Incorrect two-body acceleration");
    check_collisions(tree, particles);

    // More than one bucket's capacity, all at the same position.
    particles.assign(16, PointMass{1, {2, 2, 2}, {}, 1});
    tree.rebuild();
    require(tree.gravity(2, 0.5) == 0, "Coincident positions must have finite zero potential");
    for (const auto& p : particles) require(p.acceleration == Vec3{}, "Invalid coincident acceleration");
    check_collisions(tree, particles);
}

static void test_mass_aggregation() {
    std::vector<PointMass> particles{{2, {0, 0, 0}, {}, 0}};
    for (int i = 1; i <= 5; ++i) {
        particles.emplace_back(i, Vec3{99.0 + i, 0, 0}, Vec3{}, 0);
    }
    BarnesHut tree(particles);
    tree.rebuild();
    // The root contains the target and must open even for a large theta. The five
    // distant sources form a single cell with total mass 15 and center 1540/15.
    tree.gravity(1, 10);
    const double center = 1540.0 / 15.0;
    require(near(particles[0].acceleration[0], 15 / (center * center)),
            "Distant cell must act as a mass at its center of mass");
}

static void test_against_direct() {
    std::mt19937 gen(42);
    auto sample = [&] { return static_cast<double>(gen()) / gen.max(); };
    std::vector<PointMass> particles;
    for (int i = 0; i < 512; ++i) {
        const Vec3 pos{sample(), sample(), sample()};
        const double mass = 1 + sample();
        particles.emplace_back(mass, pos, Vec3{}, 0.03 * sample());
    }
    // Two touching particles on opposite sides of the root's split plane.
    particles.emplace_back(1, Vec3{-1, 0, 0}, Vec3{}, 0.875);
    particles.emplace_back(1, Vec3{1, 0, 0}, Vec3{}, 1.125);

    auto exact = particles;
    const double exact_energy = direct_gravity(exact, 1.5);
    BarnesHut tree(particles);
    tree.rebuild();
    check_collisions(tree, particles);

    require(near(tree.gravity(1.5, 0), exact_energy), "Theta zero disagrees with exact energy");
    for (std::size_t i = 0; i < particles.size(); ++i) {
        for (int k = 0; k < 3; ++k) {
            require(near(particles[i].acceleration[k], exact[i].acceleration[k]),
                    "Theta zero disagrees with exact acceleration");
        }
    }

    const auto acceleration_error = [&] {
        double error = 0, reference = 0;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const auto delta = particles[i].acceleration - exact[i].acceleration;
            error += dot(delta, delta);
            reference += dot(exact[i].acceleration, exact[i].acceleration);
        }
        return std::sqrt(error / reference);
    };
    const double approx_energy = tree.gravity(1.5, 0.5);
    const double approx_error = acceleration_error();
    require(std::abs(approx_energy / exact_energy - 1) < 0.01, "Approximate energy error exceeds 1%");
    require(approx_error < 0.02, "Approximate acceleration error exceeds 2%");
    tree.gravity(1.5, 0.25);
    require(acceleration_error() < approx_error, "Smaller theta should improve this system's accuracy");

    particles[0].position = {20, -10, 0};
    particles[1].position = particles[0].position;
    tree.rebuild();
    check_collisions(tree, particles);
    require(near(tree.gravity(1.5, 0), direct_gravity(particles, 1.5)),
            "Rebuilt tree disagrees with exact energy after particles move");
}

int main() try {
    test_small_systems();
    test_mass_aggregation();
    test_against_direct();
    std::cout << "Barnes-Hut tests passed.\n";
} catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
}
