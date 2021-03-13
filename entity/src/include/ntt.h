//
// Created by leonard on 25.02.21.
//

#ifndef ENTITY_NTT_H
#define ENTITY_NTT_H

#include "basic_finite_fields.h"

#include <vector>
#include <concepts>
#include <vector>



template<typename T, class FF_class> requires std::derived_from<FF_class, FiniteField<T>>
class NttEngine1 {

public:
    explicit NttEngine1(T modulus, T N) : FF(modulus), pow_omega(N), pow_phi(N), pow_omega1(N), pow_phi1(N) {
        // static_assert(std::is_base_of_v<FF_class, FiniteField<T>>, "Template class argument is not a finite field subclass");
        T omega_plain = FF.ComputePrimitiveRootOfUnity(N);
        T phi_plain = FF.ComputePrimitiveRootOfUnity(2 * N);

        omega = FF.SwitchTo2N(omega_plain);
        phi = FF.SwitchTo2N(phi_plain);
        omega1 = FF.ModExp(omega, modulus - 2);
        phi1 = FF.ModExp(phi, modulus - 2);

        for (T i = 0; i < N; ++i) {
            pow_omega[i] = FF.ModExp(omega, i);
            pow_omega1[i] = FF.ModExp(omega1, i);
            pow_phi[i] = FF.ModExp(phi, i);
            pow_phi1[i] = FF.ModExp(phi1, i);
        }
    }

private:

    FF_class FF;
    T phi, phi1;
    T omega, omega1;
    std::vector<T> pow_omega;
    std::vector<T> pow_phi;
    std::vector<T> pow_omega1;
    std::vector<T> pow_phi1;
};



#endif //ENTITY_NTT_H
