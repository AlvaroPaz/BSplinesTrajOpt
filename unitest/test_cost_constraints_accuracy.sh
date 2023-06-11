#!/bin/bash
echo "Cinvestav 2023"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo "//!----------------------------------------------------------------------------------------------!//"
echo "//!----------------Testing Accuracy in Cost Function and Constraints Derivatives-----------------!//"
echo "//!----------------------------------------------------------------------------------------------!//"

cd build/examples/accuracy_tests
./accuracy_test_02 nao_inertial_python.urdf
./accuracy_test_02 HRP2.urdf
./accuracy_test_02 atlas_cpp_reduced.urdf
