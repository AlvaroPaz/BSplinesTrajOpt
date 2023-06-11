#!/bin/bash
echo "Cinvestav 2023"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo "//!----------------------------------------------------------------------!//"
echo "//!----------------Testing numeric dynamic-objects times-----------------!//"
echo "//!----------------------------------------------------------------------!//"

cd build/examples/time_tests
./example_02 nao_inertial_python.urdf
./example_02 HRP2.urdf
./example_02 atlas.urdf
