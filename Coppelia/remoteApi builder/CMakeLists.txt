cmake_minimum_required(VERSION 2.8.3)
project(remoteApiLibrary)

add_definitions(-DNON_MATLAB_PARSING)
add_definitions(-DMAX_EXT_API_CONNECTIONS=255)
add_definitions(-DDO_NOT_USE_SHARED_MEMORY)

ADD_LIBRARY( remoteApi SHARED extApi.c extApiPlatform.c )# shared_memory.c )

#target_link_libraries( remoteApi extApi.h extApiPlatform.h simConst.h extApiInternal.h )# shared_memory.h )

# commet line 10 when compiling SHARED library