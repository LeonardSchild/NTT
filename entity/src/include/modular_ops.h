//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_MODULAR_OPS_H
#define ENTITY_MODULAR_OPS_H

#include <cstdint>

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

// TODO: In the future derive everything but modulus
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

#endif //ENTITY_MODULAR_OPS_H
