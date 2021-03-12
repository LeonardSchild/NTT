cmake_minimum_required(VERSION 3.17)
project(entity CXX C)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Winline")

include(entity/src/CMakeLists.txt)
include(entity/test/CMakeLists.txt)
