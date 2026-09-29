add_executable(Compact main.cpp)

target_link_libraries(Compact PRIVATE my_objects)
set_target_properties(Compact PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

