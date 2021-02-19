//
// Created by leonard on 18.02.21.
//

#include <iostream>
#include "ModRing128.h"

int main() {


    ModRing128 ring(0,0,0,0);

    uint128_t a,b;
    wide_mul_128(444342423424234,54234242423,&a,&b);

    std::cout << uint64_t(b >> 64) << " " << uint64_t(b % (uint128_t(1) << 64)) << std::endl;
}