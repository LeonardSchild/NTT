//
// Created by leonard on 20.03.21.
//

#ifndef ENTITY_CU_NTT_ENGINE_CUH
#define ENTITY_CU_NTT_ENGINE_CUH

#include "interfaces.h"
#include "cu_finite_field.cuh"

template<typename T>
__global__ void hadamard_product(FiniteField<T> FF, T* vec_out, T* vec_in_left, T* vec_in_right) {
    uint32_t idx = blockIdx.x * blockDim.x + blockIdx.y;
    vec_out[idx] = FF.ModMul(vec_in_left[idx], vec_in_right[idx]);
}

template<typename T_base, typename T_base2x>
class CuNttEngine : public NttEngine<T_base, CuFiniteField<T_base, T_base2x>> {

public:
    using NttEngine<T_base,  CuFiniteField<T_base, T_base2x>>::FF;
    using typename NttEngine<T_base, CuFiniteField<T_base, T_base2x>>::Direction;

    CuNttEngine(T_base modulus, T_base N) : NttEngine<T_base, CuFiniteField<T_base, T_base2x>>(modulus, N),
                                            transform_buffer(N), powers_of_omega(N), powers_of_omega1(N), powers_of_phi(N), powers_of_phi1(N), reversed_idx(N) {
        PreCompute();
    }

    void Multiply(std::vector<T_base>& out, std::vector<T_base>& lhs, std::vector<T_base>& rhs) override {
        for (uint32_t i = 0; i < this->N; ++i) {
            out[i] = FF.ModMul(lhs[i], rhs[i]);
        }
    }

    void Multiply(T_base* out, T_base* vec_in_left, T_base*  vec_in_right) {
        auto N = this->N;
        auto X_dim = N > 1024 ? N / 1024 : 1;
        auto Y_dim = N > 1024 ? 1024 : N;

        hadamard_product<<<X_dim, Y_dim>>>(out, vec_in_left, vec_in_right);
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

        cudaMalloc(&device_powers_of_omega, sizeof(T_base) * N);
        cudaMalloc(&device_powers_of_omega1, sizeof(T_base) * N);
        cudaMalloc(&device_powers_of_phi, sizeof(T_base) * N);
        cudaMalloc(&device_powers_of_phi1, sizeof(T_base) * N);
        cudaMalloc(&device_reversed_idx, sizeof(T_base) * N);

        cudaMemcpy(device_powers_of_omega, powers_of_omega.data(), sizeof(T_base) * N, cudaMemcpyHostToDevice);
        cudaMemcpy(device_powers_of_omega1, powers_of_omega1.data(), sizeof(T_base) * N, cudaMemcpyHostToDevice);
        cudaMemcpy(device_powers_of_phi, powers_of_phi.data(), sizeof(T_base) * N, cudaMemcpyHostToDevice);
        cudaMemcpy(device_powers_of_phi1, powers_of_phi1.data(), sizeof(T_base) * N, cudaMemcpyHostToDevice);
        cudaMemcpy(device_reversed_idx, reversed_idx.data(), sizeof(T_base) * N, cudaMemcpyHostToDevice);
    }



private:

    T_base *device_powers_of_omega;
    T_base *device_powers_of_omega1;
    T_base *device_powers_of_phi;
    T_base *device_powers_of_phi1;
    T_base *device_reversed_idx;

    std::vector<T_base> powers_of_omega;
    std::vector<T_base> powers_of_omega1;
    std::vector<T_base> powers_of_phi;
    std::vector<T_base> powers_of_phi1;
    std::vector<T_base> reversed_idx;

};

#endif //ENTITY_CU_NTT_ENGINE_CUH
