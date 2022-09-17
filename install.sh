#!/bin/bash  
echo "Cinvestav 2022"


echo "Cloning Eigen3 repository"
git clone --depth 1 --branch before-git-migration https://gitlab.com/libeigen/eigen.git
mv eigen eigen3

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib

echo "Building BSpline Trajectory Optimization Library"
mkdir build
cd build
cmake -Wno-dev .. -DCMAKE_BUILD_TYPE=RELEASE
make

