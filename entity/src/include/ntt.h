//
// Created by leonard on 25.02.21.
//

#ifndef ENTITY_NTT_H
#define ENTITY_NTT_H

#include "basic_finite_fields.h"
#include "interfaces.h"
#include "utils.h"

#include <vector>
/**
 * This subclass implements the functionality of a basic NttEngine. It follows the algorithm(s) described in
 * https://eprint.iacr.org/2014/646.pdf, but has not been optimised.
 * @tparam T_base Base type for the computations within this class, and the in-out type for the underlying BasicFiniteField
 * @tparam T_base2x Internal type for the BasicFiniteField
 */
template<typename T_base, typename T_base2x>
class BasicNttEngine : NttEngine<T_base, BasicFiniteField<T_base, T_base2x>> {

    using NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>::FF;
    using typename NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>::Direction;

    BasicNttEngine(T_base modulus, T_base N) : NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>(modulus, N), transform_buffer(N) {
        this->N = N;
        PreCompute();
    }

    void Forward(std::vector<T_base>& in_out) override {

        for(uint32_t i = 0; i < this->N; i++)
            in_out[i] = FF.ModMul(in_out[i], FF.ModExp(this->phi[0], i));

        ProtoTransform(in_out, Direction::Forward);
    }

    void Backward(std::vector<T_base>& in_out) override {

        ProtoTransform(in_out, Direction::Backward);

        for (uint32_t i = 0; i < this->N; i++)
            in_out[i] = FF.ModMul(in_out[i], FF.ModExp(this->phi[1], i));
    }

    void Multiply(std::vector<T_base>& out, std::vector<T_base>& lhs, std::vector<T_base>& rhs) override {
        for (uint32_t i = 0; i < this->N; ++i) {
            out[i] = FF.ModMul(lhs[i], rhs[i]);
        }
    }

    void TransformAndMultiply(std::vector<T_base>& out, std::vector<T_base>& lhs, std::vector<T_base>& rhs) override {
        Forward(lhs);
        Forward(rhs);
        Multiply(out, lhs, rhs);
    }

protected:


    void PreCompute() override {

        T_base N = this->N;
        this->logN = log2(N);
        T_base omega_t = FF.ComputePrimitiveRootOfUnity(N);
        T_base phi_t = FF.ComputePrimitiveRootOfUnity(2 * N);
        T_base N_t = FF.SwitchTo2N(N);
        T_base modulus = FF.GetModulus();

        this->omega[0] = FF.SwitchTo2N(omega_t);
        // Use Fermat's small theorem
        this->omega[1] = FF.ModExp(this->omega[0], modulus - 2);
        this->phi[0] = FF.SwitchTo2N(phi_t);
        // Use Fermat's small theorem
        this->phi[1] = FF.ModExp(this->phi[0], modulus - 2);
        this->invN = FF.ModExp(N_t, modulus - 2);

    }

    void ProtoTransform(std::vector<T_base>& in_out, Direction d) override {

        // For this implementation we assume that the index fits into a 32-bit unsigned integer

        for (uint32_t i = 0; i < this->N[0] / 2; i++)
            transform_buffer[i] = in_out[bit_reverse(i, this->logN)];

        auto twiddle_factor = this->omega[d];

        for (uint32_t i = 0; i < this->logN; i++) {
            for(uint32_t j = 0; j < this->N / 2; j++) {
                uint32_t pij = (j >> (this->logN - i)) << (this->logN - i);
                auto o = FF.ModExp(twiddle_factor, pij);
                auto a = transform_buffer[2 * j];
                auto b = FF.ModMul(transform_buffer[2 * j + 1], o);
                in_out[j] = FF.ModAdd(a,b);
                in_out[j + this->N/2] = FF.ModSub(a,b);
            }
            in_out.swap(transform_buffer);
        }

        if (this->logN & 1)
            in_out.swap(transform_buffer);

        if (d == Direction::BACKWARD) {
            for (uint32_t i = 0; i < this->N; i++) {
                in_out[i] = FF.ModMul(in_out[i], this->invN);
            }
        }
    }

private:

    std::vector<T_base> transform_buffer;

};

#endif //ENTITY_NTT_H
