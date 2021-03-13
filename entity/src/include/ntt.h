//
// Created by leonard on 25.02.21.
//

#ifndef ENTITY_NTT_H
#define ENTITY_NTT_H

#include "basic_finite_fields.h"
#include "interfaces.h"

#include <vector>

template<typename T_base, typename T_base2x, bool precompute_factors>
class BasicNttEngine : NttEngine<T_base, BasicFiniteField<T_base, T_base2x>> {
    using NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>::FF;

    BasicNttEngine(T_base modulus, T_base N) : NttEngine<T_base, BasicFiniteField<T_base, T_base2x>>(modulus, N) {
        this->N[0] = N;
        PreCompute();
    }

protected:

    void PreCompute() override {

        T_base N = this->N;
        T_base omega_t = FF.ComputePrimitiveRootOfUnity(this->N);
        T_base phi_t = FF.ComputePrimitiveRootOfUnity(2 * this->N);
        T_base N_t = FF.SwitchTo2N(N);
        T_base modulus = FF.GetModulus();

        this->omega[0] = FF.SwitchTo2N(omega_t);
        this->omega[1] = FF.ModExp(this->omega[0], modulus - 2);
        this->phi[0] = FF.SwitchTo2N(phi_t);
        this->phi[1] = FF.ModExp(this->phi[0], modulus - 2);
        this->N[1] = FF.ModExp(N_t, modulus - 2);

        if constexpr (precompute_factors) {

            powers_of_omega.reserve(N);
            powers_of_omega1.reserve(N);
            powers_of_phi.reserve(N);
            powers_of_phi1.reserve(N);

            for (T_base i = 0; i < N; i++) {
                powers_of_omega[i] = FF.ModExp(this->omega[0], i);
                powers_of_omega1[i] = FF.ModExp(this->omega[1], i);
                powers_of_phi[i] = FF.ModExp(this->phi[0], i);
                powers_of_phi1[i] = FF.ModExp(this->phi1[1], i);
            }
        }
    }

private:

    std::vector<T_base> powers_of_phi;
    std::vector<T_base> powers_of_phi1;
    std::vector<T_base> powers_of_omega;
    std::vector<T_base> powers_of_omega1;

};

#endif //ENTITY_NTT_H
