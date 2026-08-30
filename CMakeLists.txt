cmake_minimum_required(VERSION 3.23)

project(Compact)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)
add_subdirectory(results/LWE)
add_subdirectory(results/nonUniformFD)
