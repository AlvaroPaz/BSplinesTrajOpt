#!/bin/bash
echo "Cinvestav 2022"

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

