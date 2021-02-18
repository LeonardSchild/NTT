//
// Created by leonard on 18.02.21.
//

#include <iostream>
#include "src/include/modular_ops.h"


int main() {

    ModRing32 test(134215681, 130021375, 0, 0);

    uint32_t x = test.SwitchTo2N(12345);
    uint32_t y = test.SwitchTo2N(56768);
    uint32_t z = test.ModAdd(x,y);
    uint32_t b = test.ModMul(x,y);
    uint32_t a = test.SwitchFrom2N(b);
    std::cout << "Running tests... " << a << " " << b << std::endl;

}