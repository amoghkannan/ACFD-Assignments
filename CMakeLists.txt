cmake_minimum_required(VERSION 3.23)

project(Compact)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(src)

add_executable(Compact)

#target_compile_options(my_objects PRIVATE
#    -g -O0 -Wall
#)

target_link_libraries(Compact my_objects)
set_target_properties(Compact PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

