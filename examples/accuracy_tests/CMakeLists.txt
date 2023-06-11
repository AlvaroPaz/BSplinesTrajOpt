message("BUILD_EXAMPLES is active")

add_executable (accuracy_test_01 accuracy_test_01.cc)
target_link_libraries (accuracy_test_01 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable (accuracy_test_02 accuracy_test_02.cc)
target_link_libraries (accuracy_test_02 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable (accuracy_test_03 accuracy_test_03.cc)
target_link_libraries (accuracy_test_03 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

