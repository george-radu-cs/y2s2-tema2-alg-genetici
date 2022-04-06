#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

#include <random>

namespace RandomGenerator {
    inline std::mt19937_64 mt{std::random_device{}()};
    inline std::uniform_int_distribution<> oneOrZeroDistribution{0, 1};

    inline int getOneOrZero() {
        return oneOrZeroDistribution(mt);
    }

    inline long double getRandomRealNumber(double lower_bound, double upper_bound) {
        std::uniform_real_distribution<> dist{lower_bound, upper_bound};
        return dist(mt);
    }

    inline int getRandomInt(int lower_bound, int upper_bound) {
        std::uniform_int_distribution<> dist{lower_bound, upper_bound};
        return dist(mt);
    }

    inline long double getRandomNumberBetweenZeroAndOne() {
        return getRandomRealNumber(0, 1);
    }
}

#endif // RANDOM_GENERATOR_H
