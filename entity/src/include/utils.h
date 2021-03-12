//
// Created by leonard on 22.02.21.
//

#ifndef ENTITY_UTILS_H
#define ENTITY_UTILS_H

inline __uint128_t make_2u64u128(uint64_t a, uint64_t b) {
    return __uint128_t(a) << 64 | b;
}

inline __uint128_t make_4u32u128(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    return ((__uint128_t(a) << 32) | b) << 64 | ((uint128_t(c) << 32) | d);
}

void print_u128(uint128_t &a) {
    uint64_t hi = a >> 64;
    uint64_t lo = (a << 64) >> 64;
    std::cout << "HI: " << hi << " LO: " << lo << std::endl;
}

#endif //ENTITY_UTILS_H
