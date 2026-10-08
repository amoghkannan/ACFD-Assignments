add_executable(Krylov main.cpp)

target_link_libraries(Krylov PRIVATE my_objects)
set_target_properties(Krylov PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

