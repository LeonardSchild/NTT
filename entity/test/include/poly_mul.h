//
// Created by leonard on 16.03.21.
//

#ifndef ENTITY_POLY_MUL_H
#define ENTITY_POLY_MUL_H

#include <vector>
#include <cassert>
#include "src/include/interfaces.h"

template<typename T>
void polymul(FiniteField<T>& FF, std::vector<T>& out, std::vector<T>& lhs, std::vector<T>& rhs) {
    assert(lhs.size() == rhs.size());
    assert(out.size() >= 2 * lhs.size());

    for (uint32_t i = 0; i < 2 * lhs.size(); i++) {
        T temp = 0;
        for (uint32_t j = 0; j <= i; ++j) {
            if ((j >= lhs.size()) || ( (i-j) >= rhs.size()))
                continue;

            temp = FF.ModAdd(temp, FF.ModMul(lhs[j], rhs[i-j]));
        }
        out[i] = temp;
    }
}

#endif //ENTITY_POLY_MUL_H
