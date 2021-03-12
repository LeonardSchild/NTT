//
// Created by leonard on 17.02.21.
//

#ifndef ENTITY_MODULAR_OPS_H
#define ENTITY_MODULAR_OPS_H

#include <stdexcept>
#include <cstdint>
#include <climits>

/** Abstract template class to model a finite field
 * @tparam T variable type for computations
 */
template<typename T>
class FiniteField {

public:
    /**
     *  Constructor
     * @param modulus modulus of finite field
     */
    explicit FiniteField(T modulus) : modulus(modulus), bitwidth(sizeof(T) * CHAR_BIT) {};

    /**
     * Destructor
     */
    virtual ~FiniteField() = default;

    /**
     * Maps elements of the finite field to their Montgomery for modulus 2^\p bitwidth
     * @param input Element to map
     * @return \p input in Montgomery form
     */
    virtual T SwitchTo2N(T input) = 0;

    /**
     * Maps element from their Montgomery form to their representation modulo \p modulus
     * @param input Element in Montgomery form
     * @return \p input modulo \p modulus
     */
    virtual T SwitchFrom2N(T input) = 0;

    /**
     * Multiplies two elements of the finite field given in Montgomery form.
     * @param a An element from the finite field in Montgomery form
     * @param b An element from the finite field in Montgomery form
     * @return The product of \p a and \p b in Montgomery form
     */
    virtual T ModMul(T a, T b) = 0;

    /**
     * Adds two elements of the finite field given in Montgomery form
     * @param a An element from the finite field in Montgomery form
     * @param b An element from the finite field in Montgomery form
     * @return The sum of \p a and \p b
     */
    virtual T ModAdd(T a, T b) = 0;

    /**
     * Subtracts two elements of the finite field given in Montgomery form
     * @param a An element from the finite field in Montgomery form
     * @param b An element from the finite field in Montgomery form
     * @return The difference of \p a and \p b
     */
    virtual T ModSub(T a, T b) = 0;

    /**
     * Raises a element of the finite field given in Montgomery Form to a power.
     * @param a An element from the finite field in Montgomery form
     * @param e The power to which \p a should be raised
     * @return \p a raised to the \p e-th power
     */
    virtual T ModExp(T a, T e) = 0;

    /**
     * Computes a primitive root of unity given the order
     * @param order Order of the primitive root. We assume it divides \p modulus - 1
     * @return A primitive root of unity of given order
     */
    virtual T ComputePrimitiveRootOfUnity(T order) = 0;

protected:

    /**
     * Precomputes parameters based on the modulus
     */
    virtual void PreCompute() = 0;

    //! Modulus of finite field
    T modulus;
    //! Negated inverse of \modulus modulo 2^\bitwidth
    T modulus_inverse;
    //! Number of bits used by the type T
    T bitwidth;

};

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
        return a + b;
    }

    inline T_base ModSub(T_base a, T_base b) override {
        return a - b;
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

        if (group_order % order != 0) {
            throw std::invalid_argument("Provided order does not divide Q-1. No element can be found");
        }

        T_base2x result = 1;
        T_base2x init = 2;
        T_base2x exp = group_order / order;

        while (exp != 0) {
            if (exp & 1)
                result = (init * result) % this->modulus;
            exp >>= 1;
            init = (init * init) % this->modulus;
        }

        return result;
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
    T_base pow2N;
};

using FiniteField32 = BasicFiniteField<uint32_t, uint64_t>;
#ifdef __SIZEOF_INT128__
using FiniteField64 = BasicFiniteField<uint64_t, __uint128_t>;
#endif

#endif //ENTITY_MODULAR_OPS_H
