
add_executable(com_example_01 main_com_01.cpp)
target_link_libraries(com_example_01 ipopt_interface_nao_iner_com_01 ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(com_example_02 main_com_02.cpp)
target_link_libraries(com_example_02 ipopt_interface_nao_iner_com ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})
