add_executable(IBM main.cpp)

target_link_libraries(IBM PRIVATE my_objects)
set_target_properties(IBM PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

