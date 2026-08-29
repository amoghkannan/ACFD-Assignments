add_library(my_objects OBJECT
logger.cpp
utils.cpp
mesh.cpp
vars.cpp
scheme.cpp
integrator.cpp
solver.cpp
LWESolver.cpp
main.cpp
)

target_include_directories(my_objects PRIVATE ${CMAKE_SOURCE_DIR}/include)  
