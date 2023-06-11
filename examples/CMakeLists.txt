message("BUILD_EXAMPLES is active")

#add_subdirectory(examples)

add_subdirectory(examples_com)

add_subdirectory(examples_mu)

add_subdirectory(time_tests)

add_subdirectory(accuracy_tests)

add_subdirectory(time_tests/com_mu_times/examples_com_analytic)

add_subdirectory(time_tests/com_mu_times/examples_mu_analytic)

add_subdirectory(time_tests/com_mu_times/examples_com_bfgs)

add_subdirectory(time_tests/com_mu_times/examples_mu_bfgs)
