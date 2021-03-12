//
// Created by leonard on 19.02.21.
//

#include "modular_ops.h"
#include <iostream>

//TODO finish this

using uint128_t = __uint128_t;

#define MASK_LO_64 ((uint128_t(1) << 64) - 1)

__attribute__((noinline)) void wide_mul_128(uint128_t a, uint128_t b, uint128_t *hi_ab, uint128_t *lo_ab) {

    uint64_t a_lo = a & MASK_LO_64;
    uint64_t a_hi = a >> 64;
    uint64_t b_lo = b & MASK_LO_64;
    uint64_t b_hi = b >> 64;

    uint128_t ab_hihi = uint128_t(a_hi) * b_hi;
    uint128_t ab_lolo = uint128_t(a_lo) * b_lo;
    uint128_t ab_hilo = uint128_t(a_hi) * b_lo;
    uint128_t ab_lohi = uint128_t(a_lo) * b_hi;

    uint128_t ab_lo = ab_lolo;
    uint128_t ab_hi = ab_hihi;

    uint128_t ovf_mid = (ab_hilo >= (0 - ab_lohi)) & (ab_lohi > 0);
    uint128_t ab_mid = ab_hilo + ab_lohi;

    uint128_t ab_mid_lo = ab_mid << 64;
    uint128_t ab_mid_hi = ab_mid >> 64;

    uint128_t ovf_midlo = (ab_lo >= (0 - ab_mid_lo)) & (ab_mid_lo > 0);
    ab_lo += ab_mid_lo;

    ab_hi += (ovf_mid << 64) + ab_mid_hi + ovf_midlo;

    std::cout << uint64_t(ab_mid_hi) << " " << uint64_t(ovf_midlo) << std::endl;

    *hi_ab = ab_hi;
    *lo_ab = ab_lo;

}

class ModRing128 : FiniteField<uint128_t> {

public:
    ModRing128(uint128_t modulus, uint128_t modulus_inverse, uint128_t phi, uint32_t omega) : FiniteField<uint128_t>(
            modulus) {

        if ((modulus >> 127) == 1) {
            // :^)
            pow128modQ = 0 - modulus;
        } else {
            uint128_t tmp = (uint128_t(1) << 127) % modulus;
            pow128modQ = (2 * tmp) % modulus;
        }
    };

    inline uint128_t SwitchTo2N(uint128_t in) override {
        return (uint128_t(in) * uint128_t(pow128modQ)) % modulus;
    }

    inline uint128_t SwitchFrom2N(uint128_t in) override {
        return ModMul(1, in);
    }

    inline uint128_t ModMul(uint128_t a, uint128_t b) override {
        uint128_t x = uint128_t(a) * uint128_t(b);
        uint64_t s = uint64_t(x % ((uint128_t(1) << 64) - 1)) * modulus_inverse;
        uint64_t t = (x + uint128_t(s) * uint128_t(modulus)) >> 64;
        return t < modulus ? t : t - modulus;
    }

    inline uint128_t ModAdd(uint128_t a, uint128_t b) override {
        return a + b;
    }

    inline uint128_t ModSub(uint128_t a, uint128_t b) override {
        return a - b;
    }

private:
    uint64_t pow128modQ = 0;
};
