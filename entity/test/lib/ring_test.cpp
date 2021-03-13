//
// Created by leonard on 19.02.21.
//

#include "src/include/entity.h"

#include "test/include/ring_test.h"
#include "test/include/default_random_engine.h"

#include <iostream>

bool test_ring_32() {

    uint32_t modulus = 134215681;
    uint32_t order = 257;

    FiniteField32 test_ring(modulus);
    DefaultRandomEngine32 engine(0, modulus, 3.2);

    uint32_t a = engine.SampleUniformInt();
    uint32_t b = engine.SampleUniformInt();

    uint32_t expected = (uint64_t(a) * uint64_t(b)) % modulus;
    uint32_t am = test_ring.SwitchTo2N(a);
    uint32_t bm = test_ring.SwitchTo2N(b);


    uint32_t mm = test_ring.ModMul(am, bm);
    uint32_t res = test_ring.SwitchFrom2N(mm);

    uint32_t root = test_ring.ComputePrimitiveRootOfUnity(order);
    uint32_t mg_root = test_ring.SwitchTo2N(root);
    uint32_t mg_pow = test_ring.ModExp(mg_root, order);
    uint32_t pow = test_ring.SwitchFrom2N(mg_pow);

    if (expected != res) {
        std::cerr << "Failed to multiply in \\ZZ_" << modulus << std::endl;
        std::cerr << "A = " << a << ", B = " << b << ". Got " << res << ", expected " << expected << std::endl;
        return false;
    }

    if (pow != 1) {
        std::cerr << "Failed to compute root of unity in \\ZZ_" << modulus << std::endl;
        std::cerr << "r = " << root << ", r^" << order << " = " << pow << std::endl;
        return false;
    }

    return true;
}

bool test_ring_64() {

    uint64_t modulus = 17179869209;
    uint64_t order = 83;

    FiniteField64 test_ring(modulus);
    DefaultRandomEngine64 engine(0, modulus, 3.2);

    uint64_t a = engine.SampleUniformInt();
    uint64_t b = engine.SampleUniformInt();

    uint64_t expected = (uint64_t(a) * uint64_t(b)) % modulus;
    uint64_t am = test_ring.SwitchTo2N(a);
    uint64_t bm = test_ring.SwitchTo2N(b);


    uint64_t mm = test_ring.ModMul(am, bm);
    uint64_t res = test_ring.SwitchFrom2N(mm);
    uint64_t root = test_ring.ComputePrimitiveRootOfUnity(order);
    uint64_t mg_root = test_ring.SwitchTo2N(root);
    uint64_t mg_pow = test_ring.ModExp(mg_root, order);
    uint64_t pow = test_ring.SwitchFrom2N(mg_pow);

    if (expected != res) {
        std::cerr << "Failed to multiply in \\ZZ_" << modulus << std::endl;
        std::cerr << "A = " << a << ", B = " << b << ". Got " << res << ", expected " << expected << std::endl;
        return false;
    }

    if (pow != 1) {
        std::cerr << "Failed to compute root of unity in \\ZZ_" << modulus << std::endl;
        std::cerr << "r = " << root << ", r^" << order << " = " << pow << std::endl;
        return false;
    }

    return true;
}

bool exec_ring_tests() {

    bool pass = test_ring_32();
    pass = pass && test_ring_64();

    if (pass) {
        std::cerr << "All tests passed" << std::endl;
    }

    return pass;
}