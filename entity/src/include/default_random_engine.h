//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_DEFAULT_RANDOM_ENGINE_H
#define ENTITY_DEFAULT_RANDOM_ENGINE_H

#include "random_engine.h"

class DefaultRandomEngine : RandomEngine<uint32_t> {

public:

    DefaultRandomEngine(uint32_t a, uint32_t b, double sigma) : RandomEngine<uint32_t>(a,b,sigma,0.f),
            generator(rd()), idistr(a,b), ndistr(0.f, sigma) {
        length = upper_bound - lower_bound;
    }

    uint32_t SampleUniform();

    uint32_t SampleNormal();

private:

    uint32_t length;
    std::random_device rd;
    std::mt19937 generator;
    std::uniform_int_distribution<uint32_t> idistr;
    std::normal_distribution<double> ndistr;
};

#endif //ENTITY_DEFAULT_RANDOM_ENGINE_H
