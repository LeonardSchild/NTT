#include <cstdint>
#include <iostream>
#include <chrono>

#include "ntt.hpp"

#define MOD 134215681
#define POLY_DIM 1024
#define usint uint64_t

int main() {
	
	auto ntt = NTT_worker<usint, MOD, POLY_DIM>();
	
	auto p1 = RingPolynomial<usint, POLY_DIM>();
	auto p2 = RingPolynomial<usint, POLY_DIM>();
	auto pr = RingPolynomial<usint, POLY_DIM>();

	for(int i = 0; i < POLY_DIM; i++) {
		p1.coefficients[i] = i;
		p2.coefficients[i] = i;
	}

	auto start = std::chrono::high_resolution_clock::now();

	ntt.multiply(pr, p1, p2);
	
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
	
	std::cout << "Product took " << elapsed.count() << "us." << std::endl;
}
