//
// Created by leonard on 18.02.21.
//

#include <iostream>
#include <chrono>

#include "test/include/ring_test.h"
#include "test/include/ntt_test.h"

int main() {
    bool res = exec_ring_tests();
    res = res && exec_ntt_tests();




    if (res)
        std::cerr << "All tests passed" << std::endl;

}