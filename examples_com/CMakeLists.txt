message("BUILD_COM_EXAMPLES is active")

#add_executable(com_example_01 main_com_01.cpp)
#target_link_libraries(com_example_01 ipopt_interface_nao_iner_com_01 ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(com_example_02 main_com_02.cpp)
target_link_libraries(com_example_02 ipopt_interface_nao_iner_com ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

#add_executable(com_example_03 main_com_03.cpp)
#target_link_libraries(com_example_03 ipopt_interface_nao_iner_com ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

#add_executable(com_example_04 main_com_04.cpp)
#target_link_libraries(com_example_04 ipopt_interface_nao_iner_com ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

#add_executable(com_example_05 main_com_05.cpp)
#target_link_libraries(com_example_05 ipopt_interface_nao_iner_com ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

#add_executable(com_example_06 main_com_06.cpp)
#target_link_libraries(com_example_06 ipopt_interface_nao_iner_com ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})
