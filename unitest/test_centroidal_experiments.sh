#!/bin/bash
echo "Cinvestav 2023"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo "//!----------------------------------------------------------------------!//"
echo "//!---------------Running analytic Centroidal experiments----------------!//"
echo "//!----------------------------------------------------------------------!//"

cd build/examples/time_tests/com_mu_times/examples_mu_analytic
./mu_analytic_02
./mu_analytic_03
./mu_analytic_04

echo "//!----------------------------------------------------------------------!//"
echo "//!-----------------Running BFGS Centroidal experiments------------------!//"
echo "//!----------------------------------------------------------------------!//"

cd ../examples_mu_bfgs
./mu_bfgs_02
./mu_bfgs_03
./mu_bfgs_04







