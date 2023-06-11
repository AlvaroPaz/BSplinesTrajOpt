#!/bin/bash
echo "Cinvestav 2023"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo "//!----------------------------------------------------------------------!//"
echo "//!------------------Running analytic COM experiments--------------------!//"
echo "//!----------------------------------------------------------------------!//"

cd build/examples/time_tests/com_mu_times/examples_com_analytic
./com_analytic_02
./com_analytic_03
./com_analytic_04
./com_analytic_05
./com_analytic_06

echo "//!----------------------------------------------------------------------!//"
echo "//!---------------------Running BFGS COM experiments---------------------!//"
echo "//!----------------------------------------------------------------------!//"

cd ../examples_com_bfgs
./com_bfgs_02
./com_bfgs_03
./com_bfgs_04
./com_bfgs_05
./com_bfgs_06






