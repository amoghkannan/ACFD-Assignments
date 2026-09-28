add_executable(Multigrid main.cpp)

target_link_libraries(Multigrid PRIVATE my_objects)
set_target_properties(Multigrid PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

