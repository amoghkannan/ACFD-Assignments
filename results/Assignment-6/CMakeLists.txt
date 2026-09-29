add_executable(icoNS main.cpp)

target_link_libraries(icoNS PRIVATE my_objects)
set_target_properties(icoNS PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

