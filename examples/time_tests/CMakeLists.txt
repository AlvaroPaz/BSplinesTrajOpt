message("BUILD_TIME_TESTS is active")

add_executable (example_01 example_01.cc)
target_link_libraries (example_01 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable (example_02 example_02.cc)
target_link_libraries (example_02 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})
