//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_MODULAR_OPS_H
#define ENTITY_MODULAR_OPS_H

#include <cstdint>

using uint128_t = __uint128_t;

template<typename T>
class ModRing {

public:
    ModRing(T modulus, T modulus_inverse, T phi, T omega) : modulus(modulus), modulus_inverse(modulus_inverse), phi(phi), omega(omega) {};
    virtual ~ModRing() = default;
    virtual T SwitchTo2N(T in) = 0;
    virtual T SwitchFrom2N(T in) = 0;
    virtual T ModMul(T a, T b) = 0;
    virtual T ModAdd(T a, T b) = 0;
    virtual T ModSub(T a, T b) = 0;

protected:

    T modulus;
    T modulus_inverse;
    T phi;
    T omega;

};

// TODO: In the future derive everything from modulus
class ModRing32 : ModRing<uint32_t> {

public:
    ModRing32(uint32_t modulus, uint32_t modulus_inverse, uint32_t phi, uint32_t omega) : ModRing<uint32_t>(modulus, modulus_inverse, phi, omega) {
        pow32modQ = uint32_t(( uint64_t(1) << 32) % modulus);
    };

    inline uint32_t SwitchTo2N(uint32_t in) override {
        return (uint64_t(in) * uint64_t(pow32modQ)) % modulus;
    }

    inline uint32_t SwitchFrom2N(uint32_t in) override {
        return ModMul(1, in);
    }

    inline uint32_t ModMul(uint32_t a, uint32_t b) override {
        uint64_t x = uint64_t(a) * uint64_t(b);
        uint32_t s = uint32_t(x & 0xffffffff) * modulus_inverse;
        uint32_t t = (x + uint64_t(s) * uint64_t(modulus)) >> 32;
        return t < modulus ? t : t - modulus;
    }

    inline uint32_t ModAdd(uint32_t a, uint32_t b) override {
        return a + b;
    }

    inline uint32_t ModSub(uint32_t a, uint32_t b) override {
        return a - b;
    }

private:
    uint32_t pow32modQ = 0;
};

#ifdef __SIZEOF_INT128__

class ModRing64 : ModRing<uint64_t> {

public:
    ModRing64(uint64_t modulus, uint64_t modulus_inverse, uint64_t phi, uint32_t omega) : ModRing<uint64_t>(modulus, modulus_inverse, phi, omega) {
        pow64modQ = uint64_t((uint128_t(1) << 64) % modulus);
    };

    inline uint64_t SwitchTo2N(uint64_t in) override {
        return (uint128_t(in) * uint128_t(pow64modQ)) % modulus;
    }

    inline uint64_t SwitchFrom2N(uint64_t in) override {
        return ModMul(1, in);
    }

    inline uint64_t ModMul(uint64_t a, uint64_t b) override {
        uint128_t x = uint128_t(a) * uint128_t(b);
        uint64_t s = uint64_t(x % ((uint128_t(1) << 64) - 1)) * modulus_inverse;
        uint64_t t = (x + uint128_t(s) * uint128_t(modulus)) >> 64;
        return t < modulus ? t : t - modulus;
    }

    inline uint64_t ModAdd(uint64_t a, uint64_t b) override {
        return a + b;
    }

    inline uint64_t ModSub(uint64_t a, uint64_t b) override {
        return a - b;
    }

private:
    uint64_t pow64modQ = 0;
};

#endif

#endif //ENTITY_MODULAR_OPS_H
