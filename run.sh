#!/bin/bash
echo "Cinvestav 2022"

echo "---------------------------------"
echo "Running dynamic-objects times"
echo "---------------------------------"

cd build/time_tests
./example_01 lwr.urdf
./example_01 nao_inertial_python.urdf
./example_01 HRP2.urdf
./example_01 atlas.urdf

