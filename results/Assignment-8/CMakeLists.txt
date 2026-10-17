add_executable(multiphase main.cpp)

target_link_libraries(multiphase PRIVATE my_objects)
set_target_properties(multiphase PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

