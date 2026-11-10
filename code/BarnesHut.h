#ifndef BARNES_HUT_H
#define BARNES_HUT_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

using Vec3 = std::array<double, 3>;

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline Vec3 operator*(double scale, const Vec3& v) {
    return {scale * v[0], scale * v[1], scale * v[2]};
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

struct PointMass {
    double mass;
    Vec3 position;
    Vec3 velocity;
    Vec3 acceleration{};
    double radius;

    PointMass(double m, const Vec3& pos, const Vec3& vel, double r)
        : mass(m), position(pos), velocity(vel), radius(r) {}

    void update(double dt) {
        for (int k = 0; k < 3; ++k) {
            velocity[k] += acceleration[k] * dt;
            position[k] += velocity[k] * dt;
        }
    }
};

inline bool check_collision(const PointMass& a, const PointMass& b) {
    return norm(a.position - b.position) <= a.radius + b.radius;
}

// Computes accelerations and returns the total gravitational potential energy.
inline double direct_gravity(std::vector<PointMass>& particles, double G) {
    double potential = 0.0;
    for (auto& p : particles) p.acceleration = {};
    for (std::size_t i = 0; i < particles.size(); ++i) {
        auto& a = particles[i];
        for (std::size_t j = i + 1; j < particles.size(); ++j) {
            auto& b = particles[j];
            const Vec3 delta = b.position - a.position;
            const double r2 = dot(delta, delta);
            if (r2 == 0.0) continue;
            const double inv_r = 1.0 / std::sqrt(r2);
            const double scale = G / r2;
            const Vec3 direction = inv_r * delta;
            for (int k = 0; k < 3; ++k) {
                a.acceleration[k] += (scale * b.mass) * direction[k];
                b.acceleration[k] -= (scale * a.mass) * direction[k];
            }
            potential -= G * a.mass * b.mass * inv_r;
        }
    }
    return potential;
}

class BarnesHut {
public:
    explicit BarnesHut(std::vector<PointMass>& particles) : particles_(particles) {}

    // Rebuild after positions change. Keep the allocated storage between steps.
    void rebuild() {
        nodes_.clear();
        order_.resize(particles_.size());
        scratch_.resize(particles_.size());
        rank_.resize(particles_.size());
        std::iota(order_.begin(), order_.end(), std::size_t{0});
        if (order_.empty()) return;

        build_node(0, order_.size(), 0);
        for (std::size_t i = 0; i < order_.size(); ++i) rank_[order_[i]] = i;
    }

    double gravity(double G, double theta) {
        double potential = 0.0;
        const double theta_squared = theta * theta;
        const std::size_t count = particles_.size();
        // The tree and particle positions are read-only during this traversal.
        // Each iteration writes its own acceleration and contributes to the energy reduction.
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(+:potential) if(count >= 1024)
#endif
        for (std::size_t i = 0; i < count; ++i) {
            Gravity result;
            accumulate_gravity(nodes_.front(), rank_[i], particles_[i].position,
                               G, theta_squared, result);
            particles_[i].acceleration = result.acceleration;
            potential += particles_[i].mass * result.potential;
        }
        // Each particle's potential includes all the other masses.
        return 0.5 * potential;
    }

    void collisions_after(std::size_t i, std::vector<std::size_t>& matches) const {
        matches.clear();
        find_collisions(nodes_.front(), i, matches);
        // Keep the direct solver's order when resolving simultaneous contacts.
        std::sort(matches.begin(), matches.end());
    }

private:
    struct Node {
        std::size_t begin, end;
        Vec3 lower, upper;
        Vec3 center_of_mass{};
        double mass = 0.0;
        double max_radius = 0.0;
        double width_squared = 0.0;
        // Child indices are packed; zero marks an unused slot (the root is index 0).
        std::array<std::size_t, 8> children{};

        bool is_leaf() const { return children[0] == 0; }
    };

    struct Gravity {
        Vec3 acceleration{};
        double potential = 0.0;
    };

    static constexpr std::size_t leaf_capacity = 4;
    static constexpr int max_depth = 48;
    std::vector<PointMass>& particles_;
    std::vector<Node> nodes_;
    std::vector<std::size_t> order_, scratch_, rank_;

    static std::size_t octant(const Vec3& pos, const Vec3& mid) {
        return (pos[0] >= mid[0] ? 1 : 0) |
               (pos[1] >= mid[1] ? 2 : 0) |
               (pos[2] >= mid[2] ? 4 : 0);
    }

    std::size_t build_node(std::size_t begin, std::size_t end, int depth) {
        const Vec3 first = particles_[order_[begin]].position;
        Node node{begin, end, first, first};
        for (std::size_t p = begin; p < end; ++p) {
            const auto& body = particles_[order_[p]];
            node.mass += body.mass;
            node.max_radius = std::max(node.max_radius, body.radius);
            for (int k = 0; k < 3; ++k) {
                node.lower[k] = std::min(node.lower[k], body.position[k]);
                node.upper[k] = std::max(node.upper[k], body.position[k]);
                node.center_of_mass[k] += body.mass * body.position[k];
            }
        }
        const Vec3 extent = node.upper - node.lower;
        const double width = *std::max_element(extent.begin(), extent.end());
        node.width_squared = width * width;
        for (double& coord : node.center_of_mass) coord /= node.mass;

        const std::size_t index = nodes_.size();
        nodes_.push_back(node);
        // Keep coincident bodies together and bound recursion for extreme clustering.
        if (end - begin <= leaf_capacity || width == 0.0 || depth == max_depth) return index;

        // Split the node's tight bounding box into eight octants, then group indices by octant.
        Vec3 mid;
        for (int k = 0; k < 3; ++k) {
            mid[k] = node.lower[k] + 0.5 * extent[k];
            if (mid[k] == node.lower[k]) mid[k] = node.upper[k];
        }
        std::array<std::size_t, 9> offsets{};
        offsets[0] = begin;
        for (std::size_t p = begin; p < end; ++p) {
            ++offsets[1 + octant(particles_[order_[p]].position, mid)];
        }
        std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
        auto next = offsets;
        for (std::size_t p = begin; p < end; ++p) {
            scratch_[next[octant(particles_[order_[p]].position, mid)]++] = order_[p];
        }
        std::copy(scratch_.begin() + begin, scratch_.begin() + end, order_.begin() + begin);

        std::size_t slot = 0;
        for (std::size_t oct = 0; oct < 8; ++oct) {
            if (offsets[oct] == offsets[oct + 1]) continue;
            const auto child = build_node(offsets[oct], offsets[oct + 1], depth + 1);
            nodes_[index].children[slot++] = child;
        }
        return index;
    }

    static void add_gravity(const Vec3& delta, double r2, double gm, Gravity& result) {
        if (r2 == 0.0) return; // Match the direct solver for coincident bodies.
        const double inv_r = 1.0 / std::sqrt(r2);
        const double scale = gm / r2;
        for (int k = 0; k < 3; ++k) result.acceleration[k] += scale * (inv_r * delta[k]);
        result.potential -= gm * inv_r;
    }

    void accumulate_gravity(const Node& node, std::size_t target_rank, const Vec3& pos,
                            double G, double theta_squared, Gravity& result) const {
        if (node.is_leaf()) {
            for (std::size_t p = node.begin; p < node.end; ++p) {
                if (p == target_rank) continue;
                const auto& source = particles_[order_[p]];
                const Vec3 delta = source.position - pos;
                add_gravity(delta, dot(delta, delta), G * source.mass, result);
            }
            return;
        }

        // Never approximate a cell containing the target: that would add self-gravity.
        if (target_rank < node.begin || target_rank >= node.end) {
            const Vec3 delta = node.center_of_mass - pos;
            const double r2 = dot(delta, delta);
            if (node.width_squared < theta_squared * r2) {
                add_gravity(delta, r2, G * node.mass, result);
                return;
            }
        }
        for (const auto child : node.children) {
            if (child == 0) break;
            accumulate_gravity(nodes_[child], target_rank, pos, G, theta_squared, result);
        }
    }

    void find_collisions(const Node& node, std::size_t i,
                         std::vector<std::size_t>& matches) const {
        const auto& target = particles_[i];
        Vec3 offset{};
        for (int k = 0; k < 3; ++k) {
            offset[k] = std::max({node.lower[k] - target.position[k],
                                  target.position[k] - node.upper[k], 0.0});
        }
        if (norm(offset) > target.radius + node.max_radius) return;

        if (node.is_leaf()) {
            for (std::size_t p = node.begin; p < node.end; ++p) {
                const auto j = order_[p];
                if (j > i && check_collision(target, particles_[j])) matches.push_back(j);
            }
            return;
        }
        for (const auto child : node.children) {
            if (child == 0) break;
            find_collisions(nodes_[child], i, matches);
        }
    }
};

#endif
