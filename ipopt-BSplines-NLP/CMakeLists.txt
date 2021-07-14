cmake_minimum_required(VERSION 2.8)
project(ipopt_example)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(IPOPT_LIBRARY_DIRS /home/alvaropaz/MY_IPOPT_DIR/build/lib)
set(IPOPT_INCLUDE_DIRS /home/alvaropaz/MY_IPOPT_DIR/build/include/coin)

message("ipopt binary path is: " ${IPOPT_LIBRARY_DIRS} )
message("ipopt include path is: " ${IPOPT_INCLUDE_DIRS} )

include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
link_directories(${IPOPT_LIBRARY_DIRS})



list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../PQP/cmake/linux-64)
find_package(PQP)
include_directories(${PQP_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../geombd/cmake/linux-64)
find_package(geombd)
include_directories(${GeoMBD_INCLUDE_DIRS})


list(APPEND CMAKE_FIND_ROOT_PATH ${PROJECT_SOURCE_DIR}/../Coppelia/cmake/linux-64)
find_package(Coppelia)
include_directories(${Coppelia_INCLUDE_DIRS})


find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})



SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt")



add_library(ipopt_interface ipopt_interface.cpp)

target_link_libraries (ipopt_interface ${GeoMBD_LIBRARIES} ${PQP_LIBRARIES} ${Coppelia_LIBRARIES} rt)



add_executable(ipopt_example_01 main_01.cpp)
target_link_libraries(ipopt_example_01 ipopt_interface)
target_link_libraries(ipopt_example_01 ipopt)
target_link_libraries(ipopt_example_01 ${GeoMBD_LIBRARIES})
target_link_libraries(ipopt_example_01 ${PQP_LIBRARIES})

add_executable(ipopt_example_02 main_02.cpp)
target_link_libraries(ipopt_example_02 ipopt_interface)
target_link_libraries(ipopt_example_02 ipopt)
target_link_libraries(ipopt_example_02 ${GeoMBD_LIBRARIES})
target_link_libraries(ipopt_example_02 ${PQP_LIBRARIES})

add_executable(ipopt_example_03 main_03.cpp)
target_link_libraries(ipopt_example_03 ipopt_interface)
target_link_libraries(ipopt_example_03 ipopt)
target_link_libraries(ipopt_example_03 ${GeoMBD_LIBRARIES})
target_link_libraries(ipopt_example_03 ${PQP_LIBRARIES})

