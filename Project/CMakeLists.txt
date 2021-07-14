cmake_minimum_required(VERSION 2.8.3)
project(Project)

set(PACKAGE_NAME "Project")
set(PACKAGE_VERSION "1.0.0")
set(PACKAGE_SO_VERSION "1.0")
set(PACKAGE_DESCRIPTION "Static-Geometric Multibody Dynamics Library")
set(PACKAGE_AUTHOR "Alvaro Paz, Gustavo Arechavaleta")
set(PACKAGE_MAINTAINER "Alvaro Paz")


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../PQP/cmake/linux-64)
find_package(PQP)
include_directories(${PQP_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../geombd/cmake/linux-64)
find_package(geombd)
include_directories(${GeoMBD_INCLUDE_DIRS})


find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


add_subdirectory(examples)

