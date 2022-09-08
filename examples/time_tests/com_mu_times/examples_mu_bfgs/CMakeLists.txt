message("BUILD_MU_EXAMPLES is active")

add_executable(mu_bfgs_01 main_mu_01.cpp)
target_link_libraries(mu_bfgs_01 ipopt_interface_nao_iner_mu_01 ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(mu_bfgs_02 main_mu_02.cpp)
target_link_libraries(mu_bfgs_02 ipopt_interface_nao_iner_mu ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(mu_bfgs_03 main_mu_03.cpp)
target_link_libraries(mu_bfgs_03 ipopt_interface_nao_iner_mu ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(mu_bfgs_04 main_mu_04.cpp)
target_link_libraries(mu_bfgs_04 ipopt_interface_nao_iner_mu ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})
