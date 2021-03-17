//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_BASIC_FINITE_FIELDS_H
#define ENTITY_BASIC_FINITE_FIELDS_H

#include <interfaces.h>

#include <stdexcept>
#include <cstdint>
#include <climits>


/**
 * Template subclass of FiniteField for primitive types
 * @tparam T_base Base type for inputs and outputs
 * @tparam T_base2x Internal type to perform computations, should be twice as big as T_base.
 */
template<typename T_base, typename T_base2x>
class BasicFiniteField : public FiniteField<T_base> {
public:

    explicit BasicFiniteField(T_base modulus) : FiniteField<T_base>(modulus) {
        PreCompute();
    }

    inline T_base SwitchTo2N(T_base in) override {
        return (T_base2x(in) * T_base2x(pow2N)) % this->modulus;
    }

    inline T_base SwitchFrom2N(T_base in) override {
        return ModMul(1, in);
    }

    inline T_base ModMul(T_base a, T_base b) override {

        T_base2x x = T_base2x(a) * T_base2x(b);
        T_base s = T_base2x(x % (T_base2x(1) << this->bitwidth)) * this->modulus_inverse;
        T_base t = (x + T_base2x(s) * T_base2x(this->modulus)) >> this->bitwidth;

        return t < this->modulus ? t : t - this->modulus;
    }

    inline T_base ModAdd(T_base a, T_base b) override {
        T_base sum = a + b;
        return sum >= this->modulus ? sum - this->modulus : sum;
    }

    inline T_base ModSub(T_base a, T_base b) override {
        return ((this->modulus) + a) - b;
    }

    T_base ModExp(T_base a, T_base e) override {
        T_base init = pow2N;

        while (e != 0) {
            if (e & 1) {
                init = ModMul(init, a);
            }
            e >>= 1;
            a = ModMul(a, a);
        }

        return init;
    }

    T_base ComputePrimitiveRootOfUnity(T_base order) override {
        // assume modulus is prime
        auto group_order = this->modulus - 1;

        if (order == 1)
            return 1;

        if (order == group_order)
            return 2;

        if (group_order % order != 0) {
            throw std::invalid_argument("Provided order does not divide Q-1. No element can be found");
        }

        T_base2x result = 1;
        T_base2x init = 2;
        while (gcd(init, group_order) != 1)
            init++;

        T_base2x exp = group_order / order;

        while (exp != 0) {
            if (exp & 1)
                result = (init * result) % this->modulus;
            exp >>= 1;
            init = (init * init) % this->modulus;
        }

        return result;
    }

    T_base GetPow2N() {
        return pow2N;
    }

protected:

    void PreCompute() override {
        pow2N = (T_base2x(1) << this->bitwidth) % this->modulus;

        // compute inverse of modulus wrt to 2^bitwidth
        // https://crypto.stackexchange.com/questions/47493/how-to-determine-the-multiplicative-inverse-modulo-64-or-other-power-of-two
        T_base2x start = 2;
        T_base2x start_mod = 1 << start;
        T_base2x mod_inverse = 0;
        for (T_base2x i = 0; i < start_mod; i++) {
            if ((this->modulus * i) % start_mod == 1) {
                mod_inverse = i;
                break;
            }
        }

        do {
            start_mod = start_mod * start_mod;
            T_base2x prod = (this->modulus * mod_inverse) % start_mod;
            T_base2x diff = prod > 2 ? (start_mod + 2) - prod : 2 - prod;
            mod_inverse = (mod_inverse * diff) % start_mod;

        } while (start_mod < (T_base2x(1) << this->bitwidth));

        this->modulus_inverse = 0 - mod_inverse;
    }

private:

    T_base gcd(T_base a, T_base b) {
        while (b != 0) {
            T_base t = b;
            b = a % b;
            a = t;
        }
        return a;
    }


    T_base pow2N;
};

#endif //ENTITY_BASIC_FINITE_FIELDS_H
