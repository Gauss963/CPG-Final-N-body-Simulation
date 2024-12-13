#include <random>

extern "C" {
    void initialize_mt(int seed);
    double generate_random();

    std::mt19937 mt_engine;

    void initialize_mt(int seed) {
        mt_engine.seed(seed);
    }

    double generate_random() {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(mt_engine);
    }
}
