cmake_minimum_required(VERSION 2.8)
project(ipopt_example)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-msse -msse2 -msse3 -O3 -pedantic -DNDEBUG")


set(IPOPT_LIBRARY_DIRS /usr/local/lib)
set(IPOPT_INCLUDE_DIRS /usr/local/include/coin-or)

message("ipopt binary path is: " ${IPOPT_LIBRARY_DIRS} )
message("ipopt include path is: " ${IPOPT_INCLUDE_DIRS} )

include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
link_directories(${IPOPT_LIBRARY_DIRS})



## ToDo: add optional for local ipopt installation
##--------------------------------------------------------------------------
#set(IPOPT_LIBRARY_DIRS /home/$ENV{USERNAME}/MY_IPOPT_DIR/build/lib)
#set(IPOPT_INCLUDE_DIRS /home/$ENV{USERNAME}/MY_IPOPT_DIR/build/include/coin)

#message("ipopt binary path is: " ${IPOPT_LIBRARY_DIRS} )
#message("ipopt include path is: " ${IPOPT_INCLUDE_DIRS} )

#include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
#link_directories(${IPOPT_LIBRARY_DIRS})
##--------------------------------------------------------------------------


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/PQP/cmake/linux-64)
find_package(PQP)
include_directories(${PQP_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/geombd/cmake/linux-64)
find_package(geombd)
include_directories(${GeoMBD_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/Coppelia/cmake/linux-64)
find_package(Coppelia)
include_directories(${Coppelia_INCLUDE_DIRS})


set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})
message(${EIGEN3_INCLUDE_DIR})


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")


add_library(ipopt_interface ipopt/ipopt_interface.cpp)
target_link_libraries(ipopt_interface ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)


add_library(ipopt_interface_nao_iner_com_01 ipopt/ipopt_interface_nao_iner_com_01)
target_link_libraries(ipopt_interface_nao_iner_com_01 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)

add_library(ipopt_interface_nao_iner_com ipopt/ipopt_interface_nao_iner_com)
target_link_libraries(ipopt_interface_nao_iner_com ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)


add_library(ipopt_interface_nao_iner_mu_01 ipopt/ipopt_interface_nao_iner_mu_01)
target_link_libraries(ipopt_interface_nao_iner_mu_01 ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)

add_library(ipopt_interface_nao_iner_mu ipopt/ipopt_interface_nao_iner_mu)
target_link_libraries(ipopt_interface_nao_iner_mu ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)

add_library(ipopt_interface_atlas_iner_com ipopt/ipopt_interface_atlas_iner_com)
target_link_libraries(ipopt_interface_atlas_iner_com ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)

add_library(ipopt_interface_atlas_iner_mu ipopt/ipopt_interface_atlas_iner_mu)
target_link_libraries(ipopt_interface_atlas_iner_mu ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)


add_subdirectory(examples)

