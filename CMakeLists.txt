cmake_minimum_required(VERSION 3.23)
project(x_spline)

set(CMAKE_CXX_STANDARD 20)
find_package(OpenMP REQUIRED)

add_executable(x_spline main.cpp)
target_link_libraries(x_spline OpenMP::OpenMP_CXX)
