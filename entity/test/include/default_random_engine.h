//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_DEFAULT_RANDOM_ENGINE_H
#define ENTITY_DEFAULT_RANDOM_ENGINE_H

#include "random_engine.h"
#include <random>

template<typename T>
class DefaultRandomEngine : RandomEngine<T> {

public:

    DefaultRandomEngine(T a, T b, double sigma) : RandomEngine<T>(a, b, sigma, 0.f),
                                                  generator(rd()), idistr(a, b), ndistr(0.f, sigma) {
        length = this->upper_bound - this->lower_bound;
    }

    T SampleUniformInt() override {
        return idistr(generator);
    }

    T SampleNormalInt() override {
        auto normal = ndistr(generator);
        auto normal_abs = fabs(normal);
        auto u_normal_abs = uint32_t(normal_abs) % length;
        if (normal < 0) {
            return this->upper_bound - u_normal_abs;
        } else {
            return this->lower_bound + u_normal_abs;
        }
    }

    ~DefaultRandomEngine() override = default;

private:

    T length;
    std::random_device rd;
    std::mt19937 generator;
    std::uniform_int_distribution<uint32_t> idistr;
    std::normal_distribution<double> ndistr;
};

using DefaultRandomEngine32 = DefaultRandomEngine<uint32_t>;
using DefaultRandomEngine64 = DefaultRandomEngine<uint64_t>;

#endif //ENTITY_DEFAULT_RANDOM_ENGINE_H
