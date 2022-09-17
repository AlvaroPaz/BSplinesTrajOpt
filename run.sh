#!/bin/bash
echo "Cinvestav 2022"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
echo $LD_LIBRARY_PATH

echo "//!----------------------------------------------------------------------!//"
echo "-------------Running analytic and numeric dynamic-objects times-------------"
echo "//!----------------------------------------------------------------------!//"

cd build/examples/time_tests
./example_01 lwr.urdf
./example_02 lwr.urdf
./example_01 nao_inertial_python.urdf
./example_02 nao_inertial_python.urdf
./example_01 HRP2.urdf
./example_02 HRP2.urdf
./example_01 atlas.urdf
./example_02 atlas.urdf

echo "//!----------------------------------------------------------------------!//"
echo "-------------Running analytic COM examples-------------"
echo "//!----------------------------------------------------------------------!//"

cd com_mu_times/examples_com_analytic
./com_analytic_02
./com_analytic_03
./com_analytic_04
./com_analytic_05
./com_analytic_06

echo "//!----------------------------------------------------------------------!//"
echo "-------------Running BFGS COM examples-------------"
echo "//!----------------------------------------------------------------------!//"

cd ../examples_com_bfgs
./com_bfgs_02
./com_bfgs_03
./com_bfgs_04
./com_bfgs_05
./com_bfgs_06

echo "//!----------------------------------------------------------------------!//"
echo "-------------Running analytic MU examples-------------"
echo "//!----------------------------------------------------------------------!//"

cd ../examples_mu_analytic
./mu_analytic_02
./mu_analytic_03
./mu_analytic_04

echo "//!----------------------------------------------------------------------!//"
echo "-------------Running BFGS MU examples-------------"
echo "//!----------------------------------------------------------------------!//"

cd ../examples_mu_bfgs
./mu_bfgs_02
./mu_bfgs_03
./mu_bfgs_04







