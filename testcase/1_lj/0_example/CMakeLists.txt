cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_STANDARD 11)

# Tell cmake where to find the feasst library
set(CMAKE_PREFIX_PATH "../../../build")
##set(CMAKE_PREFIX_PATH "$ENV{HOME}/feasst/build/")

find_package(feasst REQUIRED)

# OMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable (main test.cc)
target_link_libraries (main LINK_PUBLIC feasst)
