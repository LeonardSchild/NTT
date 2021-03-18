//
// Created by leonard on 18.03.21.
//

#ifndef ENTITY_FAST_NTT_ENGINE_H
#define ENTITY_FAST_NTT_ENGINE_H

#include "basic_finite_fields.h"
#include "interfaces.h"
#include "utils.h"


template<typename T_base, typename T_base2x>
class FastNttEngine : public NttEngine<T_base, BasicFiniteField<T_base, T_base2x>> {

public:

    using NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>::FF;
    using typename NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>::Direction;

    FastNttEngine(T_base modulus, T_base N) : NttEngine<T_base, BasicFiniteField <T_base, T_base2x>>(modulus, N),
    transform_buffer(N), powers_of_omega(N), powers_of_omega1(N), powers_of_phi(N), powers_of_phi1(N), reversed_idx(N) {
        PreCompute();
    }

    void Forward(std::vector<T_base>& in_out) override {

        for(uint32_t i = 0; i < this->N; i++)
            in_out[i] = FF.ModMul(in_out[i], powers_of_phi[i]);

        ProtoTransform(in_out, Direction::FORWARD);
    }

    void Backward(std::vector<T_base>& in_out) override {

        ProtoTransform(in_out, Direction::BACKWARD);

        for (uint32_t i = 0; i < this->N; i++)
            in_out[i] = FF.ModMul(in_out[i], powers_of_phi1[i]);
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

    FiniteField<T_base>& GetFF() override {
        return FF;
    }

protected:

    void PreCompute() override {

        T_base N = this->N;
        this->logN = mylog2(N);
        T_base phi_t = FF.ComputePrimitiveRootOfUnity(2 * N);
        T_base N_t = FF.SwitchTo2N(N);
        T_base modulus = FF.GetModulus();

        this->phi[0] = FF.SwitchTo2N(phi_t);
        // Use Fermat's small theorem
        this->phi[1] = FF.ModExp(this->phi[0], modulus - 2);

        this->omega[0] = FF.ModMul(this->phi[0], this->phi[0]);
        // Use Fermat's small theorem
        this->omega[1] = FF.ModExp(this->omega[0], modulus - 2);

        this->invN = FF.ModExp(N_t, modulus - 2);

        for (uint32_t i = 0; i < N; ++i) {
            reversed_idx[i] = bit_reverse(i, this->logN);
            powers_of_omega[i] = FF.ModExp(this->omega[0], i);
            powers_of_omega1[i] = FF.ModExp(this->omega[1], i);
            powers_of_phi[i] = FF.ModExp(this->phi[0], i);
            powers_of_phi1[i] = FF.ModExp(this->phi[1], i);
        }
    }

    void ProtoTransform(std::vector<T_base>& in_out, Direction d) override {

        // For this implementation we assume that the index fits into a 32-bit unsigned integer

        for (uint32_t i = 0; i < this->N; i++)
            transform_buffer[i] = in_out[reversed_idx[i]];

        auto twiddle_factors = d == 0 ? powers_of_omega : powers_of_omega1;

        for (uint32_t i = 0; i < this->logN; i++) {
            uint32_t mask = (1 << (this->logN - i)) - 1;
            for (uint32_t j = 0; j < this->N / 2; j++) {
                uint32_t pij = j & mask;
                T_base o = twiddle_factors[pij];
                T_base a = transform_buffer[2 * j];
                T_base b = FF.ModMul(transform_buffer[2 * j + 1], o);
                in_out[j] = FF.ModAdd(a, b);
                in_out[j + this->N / 2] = FF.ModSub(a, b);
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
    std::vector<T_base> powers_of_omega;
    std::vector<T_base> powers_of_omega1;
    std::vector<T_base> powers_of_phi;
    std::vector<T_base> powers_of_phi1;
    std::vector<T_base> reversed_idx;

};

#endif //ENTITY_FAST_NTT_ENGINE_H
