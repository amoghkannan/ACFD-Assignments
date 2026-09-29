add_executable(nonUniformFD main.cpp)

target_link_libraries(nonUniformFD PRIVATE my_objects)
set_target_properties(nonUniformFD PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

