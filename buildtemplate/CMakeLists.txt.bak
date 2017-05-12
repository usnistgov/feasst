cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 0x)
project(feasst)

# Here are some user options to customize external libraries
option(USE_XDRFILE "Use xdrfile library" OFF)
set(XDRFILE_DIR "$ENV{HOME}")
option(USE_FFTW "Use FFTW for fourier transforms" OFF)
set(FFTW_DIR "$ENV{HOME}/software/fftw-3.3.4/build")
option(USE_HDF5 "Use HDF5" OFF)
set(HDF5_USER_DIR "/usr/local/hdf5")
option(USE_GSL "Use GSL" OFF)
set(GSL_USER_DIR "$ENV{HOME}/software/gsl-2.3")
# option(USE_JSON "Use JSON" OFF)
option(USE_GTEST "Use gtest" OFF)
option(USE_CCACHE "Use ccache to speed up builds" OFF)
option(USE_OMP "Require use of OMP for parallelization" OFF)
option(USE_SWIG "Use SWIG for python interface" OFF)

# Users typically would not edit below this comment
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")

# # JSON Library
# if (USE_JSON)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_ -std=c++11")
# endif (USE_JSON)

# GSL Library
if (USE_GSL)
  find_package (GSL)
  if (GSL_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGSL_")
    include_directories(${GSL_INCLUDE_DIRS})
    set(EXTRA_LIBS ${EXTRA_LIBS} ${GSL_LIBRARIES})
  else ()
    message("GSL not found by cmake.")
    message("Attempting to search in user-defined ${GSL_USER_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGSL_")
    include_directories("${GSL_USER_DIR}/include")
    set(EXTRA_LIBS ${EXTRA_LIBS} -L${GSL_USER_DIR}/lib -lgsl -lgslcblas -lm)
  endif()
endif (USE_GSL)

# HDF5 Library
if (USE_HDF5)
	find_package(HDF5)
	if (HDF5_FOUND)
		include_directories(${HDF5_INCLUDE_DIRS})
    set(EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARY_DIRS})
  else ()
    message("HDF5 not found by cmake.")
    message("Attempting to search in user-defined ${HDF5_USER_DIR}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHDF5_")
    include_directories("${HDF5_USER_DIR}/include")
    set(EXTRA_LIBS ${EXTRA_LIBS} -L${HDF5_USER_DIR}/lib -lhdf5 -lhdf5_cpp)
	endif (HDF5_FOUND)
endif (USE_HDF5)

# FFTW Library
if (USE_FFTW)
  message("USING FFTW")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFFTW_")
  include_directories("${FFTW_DIR}/include")
  set(EXTRA_LIBS ${EXTRA_LIBS} -L${FFTW_DIR}/lib -lfftw3)
endif (USE_FFTW)

# XDR Library
if (USE_XDRFILE)
  message("USING XDRFILE")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXDRFILE_H_ -DCPLUSPLUS")
  include_directories("${XDRFILE_DIR}/include/xdrfile")
  set(EXTRA_LIBS ${EXTRA_LIBS} -L${XDRFILE_DIR}/lib -lxdrfile)
endif (USE_XDRFILE)

# ccache
if (USE_CCACHE)
  message("USING CCACHE")
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
    message("CCACHE FOUND")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  endif(CCACHE_PROGRAM)
endif(USE_CCACHE)

# OMP
find_package(OpenMP)
if(OPENMP_FOUND)
  message("OPENMP FOUND")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  
  # custom macros to detect OMP in FEASST code
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DOMP_H_")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOMP_H_")
else(OPENMP_FOUND)
  if(USE_OMP)
    message( FATAL_ERROR "OMP Requested but not found" )
    #message( SEND_ERROR "OMP Requested but not found" )
  endif()
endif()

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=gnu++0x)
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  add_definitions(-std=c++0x)
  #add_definitions(-std=c++0x -libstdc++)
endif()

file(GLOB feasst_unittest_SRC "${CMAKE_SOURCE_DIR}/../src/*_unittest.cc")
file(GLOB feasst_SRC "${CMAKE_SOURCE_DIR}/../src/*.h" "${CMAKE_SOURCE_DIR}/../src/*.cc")

# remove the unittests and main from feasst_SRC
list(REMOVE_ITEM feasst_SRC ${feasst_unittest_SRC})
file(GLOB main_SRC "${CMAKE_SOURCE_DIR}/../src/main.cc" "${CMAKE_SOURCE_DIR}/../src/ui_text.cc")
list(REMOVE_ITEM feasst_SRC ${main_SRC})

set (feasst_VERSION_MAJOR 0)
set (feasst_VERSION_MINOR 1)

include_directories("${CMAKE_SOURCE_DIR}/../src/")
if(EXISTS "${CMAKE_SOURCE_DIR}/../drivers/main.cc")
  add_executable(main "${CMAKE_SOURCE_DIR}/../drivers/main" ${feasst_SRC})
  target_link_libraries(main ${EXTRA_LIBS})
endif()
add_executable(ui_text "${CMAKE_SOURCE_DIR}/../src/ui_text.cc" ${feasst_SRC})
target_link_libraries(ui_text ${EXTRA_LIBS})
set(EXECUTABLE_OUTPUT_PATH "bin")

if (USE_GTEST)
  message("USING GTEST")
  enable_testing()
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/tmp")
  # Download and unpack googletest at configure time
  configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                   ${CMAKE_BINARY_DIR}/googletest-build)

  # The gtest/gtest_main targets carry header search path
  # dependencies automatically when using CMake 2.8.11 or
  # later. Otherwise we have to add them here ourselves.
  if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
  endif()

  # Now simply link against gtest or gtest_main as needed. Eg
  add_executable(unittest ${feasst_unittest_SRC} ${feasst_SRC})
  target_link_libraries(unittest gtest_main ${EXTRA_LIBS})
  add_test(NAME unittest COMMAND unittest)
endif(USE_GTEST)

# SWIG
if (USE_SWIG)
  FIND_PACKAGE(SWIG)
  if (SWIG_FOUND)
    message("FOUND SWIG")
    INCLUDE(${SWIG_USE_FILE})

    FIND_PACKAGE(PythonLibs)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

    #INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    SET(CMAKE_SWIG_FLAGS "")

    file(GLOB swig_SRC "${CMAKE_SOURCE_DIR}/../src/*.i")
    file(GLOB feasst_cc_SRC "${CMAKE_SOURCE_DIR}/../src/*.cc")
    list(REMOVE_ITEM feasst_cc_SRC ${feasst_unittest_SRC})
    list(REMOVE_ITEM feasst_cc_SRC ${main_SRC})

    SET_SOURCE_FILES_PROPERTIES(${swig_SRC} PROPERTIES CPLUSPLUS ON)
    #SET_SOURCE_FILES_PROPERTIES(${swig_SRC} PROPERTIES SWIG_FLAGS "-includeall")
    SWIG_ADD_MODULE(feasst python ${swig_SRC} ${feasst_cc_SRC})
    SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS})
    add_test(NAME swigtest COMMAND ../src/swigtest/swigtest.py)
  else (SWIG_FOUND)
    message("SWIG NOT FOUND")
  endif (SWIG_FOUND)
endif(USE_SWIG)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git describe --abbrev=10 --dirty --always --tags
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_definitions("-DVERSION=\"\\\"${GIT_COMMIT_HASH} ${GIT_BRANCH}\\\"\"")
add_definitions("-DFEASST_SRC_=${CMAKE_SOURCE_DIR}/../src/")
