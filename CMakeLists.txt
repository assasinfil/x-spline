cmake_minimum_required(VERSION 3.19)
project(x_spline)

set(CMAKE_CXX_STANDARD 14)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

include_directories(.)

add_executable(x_spline
        main.cpp)
