cmake_minimum_required(VERSION 3.10)
project(parallelspice)
set(CMAKE_PREFIX_PATH cspice)

#findOpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(FATAL_ERROR "Did not find OpenMP package")
endif()

# get cspice headers location
find_path(CSPICE_INCLUDE_DIR NAMES SpiceUsr.h)
if(NOT CSPICE_INCLUDE_DIR)
    message(FATAL_ERROR "Did not find SpiceUsr.h")
endif()

# get cspice library location
find_library(CSPICE_LIBRARY NAMES cspice)
if(NOT CSPICE_LIBRARY)
message(FATAL_ERROR "Did not find CSpice library")
endif()

#include cspice headers
include_directories(${CSPICE_INCLUDE_DIR})

#generate app
add_executable(parallelspice parallelspice.cpp ${HEADERS})

#link against cspice lib
target_link_libraries(parallelspice ${CSPICE_LIBRARY} OpenMP::OpenMP_CXX)
