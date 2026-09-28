add_executable(Precon main.cpp)

target_link_libraries(Precon PRIVATE my_objects)
set_target_properties(Precon PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

