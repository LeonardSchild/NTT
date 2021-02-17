//
// Created by leonard on 17.02.21.
//
#include <cstdint>
#include <random>

#include "default_random_engine.h"

uint32_t DefaultRandomEngine::SampleUniform() {
    return idistr(generator);
}

uint32_t DefaultRandomEngine::SampleNormal() {
    auto normal = ndistr(generator);
    auto normal_abs = fabs(normal);
    auto u_normal_abs = uint32_t(normal_abs) % length;
    if (normal < 0) {
        return upper_bound - u_normal_abs;
    } else {
        return lower_bound + u_normal_abs;
    }
}
