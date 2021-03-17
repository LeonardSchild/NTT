cmake_minimum_required(VERSION 3.17)
project(entity CXX C)

set(CMAKE_CXX_STANDARD 20)

set(DEBUG ON)

if (DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Winline -std=c++2a -fsanitize=address")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++2a -O3 -march=native")
endif ()


include(entity/src/CMakeLists.txt)
include(entity/test/CMakeLists.txt)
