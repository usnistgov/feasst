# If running the executable results in the error that a shared object file cannot be found,
# export LD_LIBRARY_PATH=/path/to/feasst/build/

cmake_minimum_required(VERSION 3.5)
project (FEASST_APP)
set(CMAKE_CXX_STANDARD 14)

# Tell CMake where to find the feasst library
set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR}/../build/)
#set(CMAKE_PREFIX_PATH "$ENV{HOME}/feasst/build/")
#message(${CMAKE_PREFIX_PATH})

find_package(feasst REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -g -pthread")

# set the C++14 standard
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
  endif ()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_definitions(-std=c++0x)
  endif ()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    add_definitions(-std=c++14)
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 14)
endif ()

# OMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if (NOT DEFINED FEASST_VERBOSE_LEVEL)
  set (FEASST_VERBOSE_LEVEL "3")
endif()
add_definitions("-DFEASST_VERBOSE_LEVEL_=${FEASST_VERBOSE_LEVEL}")
if (NOT DEFINED FEASST_DIR)
  set (FEASST_DIR ${CMAKE_SOURCE_DIR})
endif()
add_definitions("-DFEASST_DIR_=\"${FEASST_DIR}\"")
if (NOT DEFINED FEASST_VERSION)
  set (FEASST_VERSION "custom")
endif()
add_definitions("-DFEASST_VERSION_=\"${FEASST_VERSION}\"")

include_directories("${CMAKE_PREFIX_PATH}/include")
add_executable (tutorial tutorial.cpp)
find_package(MPI)
if (MPI_FOUND)
  target_link_libraries (tutorial LINK_PUBLIC feasstlib MPI::MPI_CXX)
else (MPI_FOUND)
  message("not found")
  target_link_libraries (tutorial LINK_PUBLIC feasstlib)
endif (MPI_FOUND)
