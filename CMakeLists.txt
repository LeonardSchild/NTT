cmake_minimum_required(VERSION 3.17)
project(entity CXX C)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Winline -std=c++2a")

include(entity/src/CMakeLists.txt)
include(entity/test/CMakeLists.txt)
