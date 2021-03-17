//
// Created by leonard on 16.03.21.
//

#include "test/include/ntt_test.h"
#include "test/include/poly_mul.h"
#include "src/include/entity.h"

bool test_basic_ntt32() {

    uint32_t modulus = 134215681;
    uint32_t poly_dim = 1 << 10;
    bool result = true;

    BasicNttEngine32 engine32(modulus, poly_dim);
    auto& FF32 = engine32.GetFF();

    std::vector<uint32_t> A(poly_dim);
    std::vector<uint32_t> B(poly_dim);
    std::vector<uint32_t> C(poly_dim);
    std::vector<uint32_t> D(2 * poly_dim);

    for (uint32_t i = 0; i < poly_dim; i++) {
        A[i] = FF32.SwitchTo2N(i);
        B[i] = FF32.SwitchTo2N(i);
    }

    polymul(FF32, D,A,B);

    for (uint32_t i = 0; i < poly_dim; ++i) {
        D[i] = FF32.ModSub(D[i], D[i + poly_dim]);
    }

    engine32.TransformAndMultiply(C,A,B);
    engine32.Backward(C);
    engine32.Backward(A);

    for (uint32_t i = 0; i < poly_dim; ++i) {
        result = result && (FF32.SwitchFrom2N(A[i]) == i);
    }

    return result;
}

bool exec_ntt_tests() {
    bool res = test_basic_ntt32();

    return res;
}