//
// Created by leonard on 18.03.21.
//

#include <chrono>
#include <vector>

#include "test/include/benchmark.h"
#include "src/include/entity.h"

template<typename T, class FF_class>
uint64_t run_benchmark(NttEngine<T, FF_class>& Ntt, uint32_t deg) {

    std::vector<T> A(deg);
    std::vector<T> B(deg);
    std::vector<T> C(deg);

    for (uint32_t i = 0; i < deg; i++) {
        A[i] = T(i);
        B[i] = T((deg + i));
    }

    auto start = std::chrono::high_resolution_clock::now();

    //Ntt.Forward(A);
    //Ntt.Forward(B);
    //Ntt.Multiply(C,A,B);
    Ntt.Backward(C);

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);

    return elapsed.count();
}

void run_benchmarks() {

    uint32_t modulus = 134215681;
    uint32_t poly_deg = 1 << 10;
    uint32_t rounds = 100;

    BasicNttEngine32 bNTT32 = BasicNttEngine32(modulus, poly_deg);
    BasicNttEngine64 bNTT64 = BasicNttEngine64(modulus, poly_deg);

    FastNttEngine32 fNTT32 = FastNttEngine32(modulus, poly_deg);
    FastNttEngine64 fNTT64 = FastNttEngine64(modulus, poly_deg);



    uint64_t b32_time = 0;
    uint64_t b64_time = 0;
    uint64_t f32_time = 0;
    uint64_t f64_time = 0;

    for(uint32_t i = 0; i < rounds; i++) {
        b32_time += run_benchmark(bNTT32, poly_deg);
        b64_time += run_benchmark(bNTT64, poly_deg);
        f32_time += run_benchmark(fNTT32, poly_deg);
        f64_time += run_benchmark(fNTT64, poly_deg);
    }

    std::cout << "[32] Basic engine took " << b32_time / rounds << "us. " << std::endl;
    std::cout << "[64] Basic engine took " << b64_time / rounds << "us. " << std::endl;
    std::cout << "[32] Fast engine took " << f32_time / rounds << "us. " << std::endl;
    std::cout << "[64] Fast engine took " << f64_time /rounds << "us. " << std::endl;

}