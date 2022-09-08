message("BUILD_MU_EXAMPLES is active")

add_executable(mu_analytic_01 main_mu_01.cpp)
target_link_libraries(mu_analytic_01 ipopt_interface_nao_iner_mu_01 ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(mu_analytic_02 main_mu_02.cpp)
target_link_libraries(mu_analytic_02 ipopt_interface_nao_iner_mu ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(mu_analytic_03 main_mu_03.cpp)
target_link_libraries(mu_analytic_03 ipopt_interface_nao_iner_mu ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})

add_executable(mu_analytic_04 main_mu_04.cpp)
target_link_libraries(mu_analytic_04 ipopt_interface_nao_iner_mu ipopt ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES})
