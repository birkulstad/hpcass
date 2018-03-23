cmake_minimum_required(VERSION 3.9)
project(hpcass)

set(CMAKE_CXX_STANDARD 11)

add_executable(hpcass main.cpp io.h io.cpp init.cpp init.h oper.cpp oper.h)

find_package(BLAS REQUIRED)

target_link_libraries(hpcass ${BLAS_LIBRARIES})

find_package(LAPACK REQUIRED)

target_link_libraries(hpcass ${LAPACK_LIBRARIES})

find_package(MPI REQUIRED)

target_link_libraries(hpcass ${MPI_LIBRARIES})