cmake_minimum_required(VERSION 3.9)
project(hpcass)

set(CMAKE_CXX_STANDARD 11)

add_executable(hpcass main.cpp)

find_package(BLAS REQUIRED)

target_link_libraries(hpcass ${BLAS_LIBRARIES})