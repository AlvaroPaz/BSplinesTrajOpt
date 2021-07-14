cmake_minimum_required(VERSION 2.8)
project(just4test)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../PQP/cmake/linux-64)
find_package(PQP)
include_directories(${PQP_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../openhrc/cmake/linux-64)
find_package(openhrc)
include_directories(${OpenHRC_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../Coppelia/cmake/linux-64)
find_package(Coppelia)
include_directories(${Coppelia_INCLUDE_DIRS})


find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})


#add_definitions(-DNON_MATLAB_PARSING)
#add_definitions(-DMAX_EXT_API_CONNECTIONS=255)
#add_definitions(-DDO_NOT_USE_SHARED_MEMORY)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")


add_executable(example main.cpp)

target_link_libraries (example ${OpenHRC_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)

