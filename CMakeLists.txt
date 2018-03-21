cmake_minimum_required(VERSION 3.9)
project(hpcass)

set(CMAKE_CXX_STANDARD 11)

add_executable(hpcass main.cpp)
add_executable(helmholtz cg_helmholtz_sb.cpp)

find_package(BLAS REQUIRED)

target_link_libraries(hpcass ${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)

target_link_libraries(hpcass ${LAPACK_LIBRARIES})