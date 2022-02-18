message("BUILD_EXAMPLES is active")

add_executable(ipopt_example_01 main_01.cpp)
target_link_libraries(ipopt_example_01 ipopt_interface)
target_link_libraries(ipopt_example_01 ipopt)
target_link_libraries(ipopt_example_01 ${GeoMBD_LIBRARIES})
target_link_libraries(ipopt_example_01 ${PQP_LIBRARIES})

add_executable(ipopt_example_02 main_02.cpp)
target_link_libraries(ipopt_example_02 ipopt_interface)
target_link_libraries(ipopt_example_02 ipopt)
target_link_libraries(ipopt_example_02 ${GeoMBD_LIBRARIES})
target_link_libraries(ipopt_example_02 ${PQP_LIBRARIES})
