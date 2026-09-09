add_library(my_objects OBJECT
logger.cpp
utils.cpp
vars.cpp
mesh.cpp
scheme.cpp
integrator.cpp
solver.cpp
LWESolver.cpp
)

target_include_directories(my_objects PUBLIC ${CMAKE_SOURCE_DIR}/include)  
