cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 0x)
project(feasst)

# Here are some user options to customize external libraries
option(USE_XDRFILE "Use xdrfile library" OFF)
set(XDRFILE_DIR "$ENV{HOME}/software/xdrfile-1.1b/build")
option(USE_FFTW "Use FFTW for fourier transforms" OFF)
set(FFTW_DIR "$ENV{HOME}/software/fftw-3.3.4/build")
option(USE_HDF5 "Use HDF5" OFF)
set(HDF5_USER_DIR "/usr/local/hdf5")
option(USE_GSL "Use GSL" OFF)
set(GSL_USER_DIR "$ENV{HOME}/software/gsl-2.3")
# option(USE_JSON "Use JSON" OFF)
option(USE_GTEST "Use gtest" OFF)
option(USE_GTESTCASES "Use testcases" OFF)
option(USE_CCACHE "Use ccache to speed up builds" OFF)
option(USE_OMP "Require use of OMP for parallelization" OFF)
option(USE_GCOV "Use GCOV for coverage testing" OFF)
option(USE_SPHINX "Use SPHINX for documentation" OFF)
option(USE_SEGERROR "Errors result in SEGFAULT for gdb backtrace" OFF)

option(USE_SWIG "Use SWIG for python interface" OFF)
# If you would like to use the python interface, you may need to specify the
# python library and include path variables below. If not specified then
# CMake attempts to find the libraries on its own.
#set(PYTHON_INCLUDE_PATH "$ENV{HOME}/anaconda2/include/python2.7")
#set(PYTHON_LIBRARIES "$ENV{HOME}/anaconda2/lib/libpython2.7.so")
set(PYTHON_INCLUDE_PATH "$ENV{HOME}/anaconda3/include/python3.6m/")
set(PYTHON_LIBRARIES "$ENV{HOME}/anaconda3/lib/libpython3.6m.so")

# Users typically would not edit below this comment
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")  #profiling

# enable dependencies
# if you're using SPHINX, you must use GTESTCASES
if (USE_SPHINX)
  set(USE_GTESTCASES ON)
endif (USE_SPHINX)
# if you're using GTESTCASES, you must use GTEST
if (USE_GTESTCASES)
  set(USE_GTEST ON)
endif (USE_GTESTCASES)

# # JSON Library
# if (USE_JSON)
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DJSON_ -std=c++11")
# endif (USE_JSON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} modules)

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

# segfault errors
if (USE_SEGERROR)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEGERROR_")
endif (USE_SEGERROR)

# SPHINX
if (USE_SPHINX)
  find_package(Doxygen REQUIRED)
  if(DOXYGEN_FOUND)
    #message("Doxygen found")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doxygen
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
  endif(DOXYGEN_FOUND)

  find_package(Sphinx REQUIRED)

#  if(NOT DEFINED SPHINX_THEME)
      set(SPHINX_THEME "classic")
#      set(SPHINX_THEME "sphinx_rtd_theme")
#  endif()

#  if(NOT DEFINED SPHINX_THEME_DIR)
#      set(SPHINX_THEME_DIR "sphinx_rtd_theme.get_html_theme_path()")
#  endif()

  # configured documentation tools and intermediate build results
  set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  # LATEX output directory
  set(SPHINX_LATEX_DIR "${CMAKE_CURRENT_BINARY_DIR}/latex")

  # FEASST version
  set(FEASST_VERSION ${GIT_COMMIT_HASH})

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in"
    "${BINARY_BUILD_DIR}/conf.py"
    @ONLY)

  add_custom_target(html ALL
    ${SPHINX_EXECUTABLE}
        -q -b html
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

  add_dependencies(html doxygen)

  add_custom_target(latex ALL
    ${SPHINX_EXECUTABLE}
        -q -b latex
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${SPHINX_LATEX_DIR}"
    COMMENT "Building LATEX documentation with Sphinx")

  add_dependencies(latex doxygen)

  add_custom_target(pdf ALL
    ${CMAKE_COMMAND} -E chdir "${SPHINX_LATEX_DIR}" make
    COMMENT "Building PDF from LATEX documentation with Sphinx")

  add_dependencies(pdf latex)

endif (USE_SPHINX)

# GCOV
if (USE_GCOV)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fprofile-arcs -ftest-coverage")
  if(CMAKE_COMPILER_IS_GNUCXX)
      include(CodeCoverage)
      setup_target_for_coverage(coverage unittest coverage)
  endif()
endif (USE_GCOV)

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

# Use the makeFactory.sh script to make feasst.i and factories.cc
execute_process(
  COMMAND ./makeFactory.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

file(GLOB feasst_unittest_SRC "${CMAKE_SOURCE_DIR}/../src/*_unittest.cc")
file(GLOB feasst_SRC "${CMAKE_SOURCE_DIR}/../src/*.h" "${CMAKE_SOURCE_DIR}/../src/*.cc")

# remove the unittests and main from feasst_SRC
list(REMOVE_ITEM feasst_SRC ${feasst_unittest_SRC})
file(GLOB main_SRC "${CMAKE_SOURCE_DIR}/../src/main.cc" "${CMAKE_SOURCE_DIR}/../src/ui_text.cc")
list(REMOVE_ITEM feasst_SRC ${main_SRC})

include_directories("${CMAKE_SOURCE_DIR}/../src/")
if(EXISTS "${CMAKE_SOURCE_DIR}/../drivers/main.cc")
  add_executable(main "${CMAKE_SOURCE_DIR}/../drivers/main" ${feasst_SRC})
  SET_TARGET_PROPERTIES(main PROPERTIES COMPILE_FLAGS -DFEASST_NAMESPACE_)
  target_link_libraries(main ${EXTRA_LIBS})
endif()
add_executable(ui_text "${CMAKE_SOURCE_DIR}/../src/ui_text.cc" ${feasst_SRC})
SET_TARGET_PROPERTIES(ui_text PROPERTIES COMPILE_FLAGS -DFEASST_NAMESPACE_)
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
  SET_TARGET_PROPERTIES(unittest PROPERTIES COMPILE_FLAGS -DFEASST_NAMESPACE_)
  target_link_libraries(unittest gtest_main ${EXTRA_LIBS})
  add_test(NAME unittest COMMAND ./bin/unittest --gtest_shuffle)

  if (USE_GTESTCASES)
    message("USING GTEST CASES")
    execute_process(
      COMMAND ./genGTestCases.sh
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
    file(GLOB feasst_testcase_SRC "${CMAKE_SOURCE_DIR}/testcase/*.cc")
    add_executable(testcase ${feasst_testcase_SRC} ${feasst_SRC})
    SET_TARGET_PROPERTIES(testcase PROPERTIES COMPILE_FLAGS -DFEASST_NAMESPACE_)
    target_link_libraries(testcase gtest_main ${EXTRA_LIBS})
    add_test(NAME testcase COMMAND ./bin/testcase --gtest_shuffle)
  endif (USE_GTESTCASES)

endif(USE_GTEST)

# SWIG
if (USE_SWIG)
  FIND_PACKAGE(SWIG)
  if (SWIG_FOUND)
    message("FOUND SWIG")
    INCLUDE(${SWIG_USE_FILE})

    # Find Python libraries. First, check if libraries are set manually
    if (PYTHON_INCLUDE_PATH)
    else (PYTHON_INCLUDE_PATH)
      FIND_PACKAGE(PythonLibs)
    endif (PYTHON_INCLUDE_PATH)
    message("FOUND PYTHON")
    message(${PYTHON_INCLUDE_PATH})
    message(${PYTHON_LIBRARIES})
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

    #INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    SET(CMAKE_SWIG_FLAGS "")

    file(GLOB feasst_cc_SRC "${CMAKE_SOURCE_DIR}/../src/*.cc")
    list(REMOVE_ITEM feasst_cc_SRC ${feasst_unittest_SRC})
    list(REMOVE_ITEM feasst_cc_SRC ${main_SRC})

    SET_SOURCE_FILES_PROPERTIES(feasst.i PROPERTIES CPLUSPLUS ON)
    #SET_SOURCE_FILES_PROPERTIES(feasst.i PROPERTIES SWIG_FLAGS "-includeall")
    SWIG_ADD_MODULE(feasst python feasst.i ${feasst_cc_SRC})
    SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS})
    add_test(NAME swigtest COMMAND ../tools/run.sh testcase/test.py)
  else (SWIG_FOUND)
    message("SWIG NOT FOUND")
  endif (SWIG_FOUND)
endif(USE_SWIG)

