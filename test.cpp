#include <cstdint>
#include <iostream>
#include <chrono>

#include "ntt.h"

#define MOD 134215681
#define POLY_DIM 1024
#define usint uint32_t

int main() {

    auto ntt = NTT_worker<usint, MOD, POLY_DIM>();

    auto p1 = RingPolynomial<usint, POLY_DIM>();
    auto p2 = RingPolynomial<usint, POLY_DIM>();
    auto pr = RingPolynomial<usint, POLY_DIM>();

    for(int i = 0; i < POLY_DIM; i++) {
        p1.coefficients[i] = 1;
        p2.coefficients[i] = 1;
    }


    uint32_t total = 0;
    uint32_t rounds = 100;
    for(int i = 0; i < rounds; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        ntt.multiply(pr.coefficients, p1.coefficients, p2.coefficients);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        total += elapsed.count();
    }

    std::cout << "Average : " << total / rounds << "us." << std::endl;
}
