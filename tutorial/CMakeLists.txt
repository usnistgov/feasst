cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)

# Tell cmake where to find the feasst library
set(CMAKE_PREFIX_PATH "$ENV{HOME}/feasst/build/")

find_package(feasst REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -g -pthread")

# set the C++11
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
  endif ()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_definitions(-std=c++0x)
  endif ()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    add_definitions(-std=c++11)
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()

# OMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories("${CMAKE_PREFIX_PATH}/include")
add_executable (tutorial tutorial.cpp)
target_link_libraries (tutorial LINK_PUBLIC feasst)
