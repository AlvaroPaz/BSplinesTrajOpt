#!/bin/bash
echo "Cinvestav 2023"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo "//!----------------------------------------------------------------------!//"
echo "//!----------------Testing analytic dynamic-objects times----------------!//"
echo "//!----------------------------------------------------------------------!//"

cd build/examples/time_tests
./example_01 nao_inertial_python.urdf
./example_01 HRP2.urdf
./example_01 atlas.urdf
