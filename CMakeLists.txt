cmake_minimum_required (VERSION 3.5)
project (FEASST)
set (CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0042 NEW) # removes warning about macosx_rpath
set(FEASST_VERSION 0.25.11)

# set path to XDR if using xtc plugin
set(XDRFILE_DIR "$ENV{HOME}/software/xdrfile-1.1b/build")

option(EV "Enable developer tools (feasst/dev/tools/depend.py and compiler debug)" OFF)
option(USE_GTEST "Use gtest" OFF)
option(USE_CCACHE "Use ccache to speed up builds" ON)
option(USE_OPENMP "Use OpenMP" ON)
option(USE_SPHINX "Use SPHINX for documentation" OFF)
#option(USE_SWIG "Use SWIG for python interface" OFF)
option(USE_GCOV "Use Coverage" OFF)
option(USE_FFTW "Use FFTW" OFF)
option(USE_PYBIND11 "Use PYBIND11" OFF)
option(USE_MPI "Use MPI" OFF)
option(USE_HEADER_CHECK "Use stand-alone header check (requires cleanup)" OFF)
set(FFTW_DIR "$ENV{HOME}/software/fftw-3.3.10/build/")
option(USE_NETCDF "Use NetCDF" OFF)
set(NETCDF_DIR "$ENV{HOME}/local")
if (NOT DEFINED FEASST_VERBOSE_LEVEL)
  set (FEASST_VERBOSE_LEVEL "3")
endif()

# If you would like to use the python interface, you may need to specify the
# python library and include path variables below. If not specified then
# CMake attempts to find the libraries on its own.
option(SET_PYTHON_PATH "Manually define PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" OFF)

if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
  add_compile_definitions(IS_X86)
endif ()

# List the plugins to compile.
if (NOT FEASST_PLUGINS)
  #set(FEASST_PLUGINS "threads")
  #set(FEASST_PLUGINS "threads;utils")
  #set(FEASST_PLUGINS "threads;utils;math")
  #set(FEASST_PLUGINS "threads;utils;math;configuration")
  #set(FEASST_PLUGINS "threads;utils;math;configuration;system")
  #set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;models;monte_carlo")
  set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;mayer;patch;example;shape;confinement;chain;flat_histogram;charge;morph;cluster;prefetch;egce;beta_expanded;aniso;gibbs;model_expanded;server")
else()
  # Note: "cmake -DEV=ON" must be enabled if FEASST_PLUGINS are changed in order to update feasst.h
  set(EV ON)
endif()

if (EV)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -g -pthread")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -g -pthread -fanalyzer")
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pg -g") #profile
  message("Enabling FEASST Developer Tools.")
else (EV)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread")
endif (EV)

# MPI
if (USE_MPI)
  set(FEASST_PLUGINS "${FEASST_PLUGINS};mpi")
endif (USE_MPI)

# print plugins
string (REPLACE ";" " " FEASST_PLUGINS_STR "${FEASST_PLUGINS}")
message(${FEASST_PLUGINS_STR})
file(WRITE "${CMAKE_BINARY_DIR}/plugins.txt" "${FEASST_PLUGINS_STR}")

# set the install location (default: in the binary build directory)"
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

# CMAKE plugins in the cmake directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/dev/cmake)

# XDR Library
list (FIND FEASST_PLUGINS "xtc" _index)
if (${_index} GREATER -1)
  message("USING XDRFILE")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DXDRFILE_H_ -DCPLUSPLUS")
  include_directories("${XDRFILE_DIR}/include/xdrfile")
  set(EXTRA_LIBS "${EXTRA_LIBS} -L${XDRFILE_DIR}/lib -lxdrfile")
endif()

# GCOV
if (USE_GCOV)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fprofile-arcs -ftest-coverage")
  include(CodeCoverage)
  setup_target_for_coverage(coverage unittest coverage --gtest_filter=-*LONG)
endif (USE_GCOV)

# FFTW
if (USE_FFTW)
#  message("USING FFTW")
#  find_package(PkgConfig REQUIRED)
#  pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
#  include_directories(PkgConfig::FFTW)
#  link_libraries(PkgConfig::FFTW)
  configure_file(${CMAKE_SOURCE_DIR}/dev/cmake/downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
  else()
    message("CMake step for findFFTW completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
  endif()

  set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
  set(CMAKE_PREFIX_PATH "${FFTW_DIR}")
  find_package(FFTW)
  if(FFTW_FOUND)
    message("FFTW FOUND")
    include_directories(${FFTW_INCLUDE_DIRS})
    link_libraries(${FFTW_LIBRARIES})
  else(FFTW_FOUND)
    message("FFTW NOT FOUND")
  endif(FFTW_FOUND)
endif (USE_FFTW)

# NetCDF
if (USE_NETCDF)
  message("USING NETCDF")
  find_package(PkgConfig REQUIRED)
  set(CMAKE_PREFIX_PATH "${NETCDF_DIR}")
  pkg_search_module(ZLIB REQUIRED zlib IMPORTED_TARGET)
  link_libraries(PkgConfig::ZLIB)
  pkg_search_module(NETCDF REQUIRED netcdf IMPORTED_TARGET)
  link_libraries(PkgConfig::NETCDF)
  pkg_search_module(NETCDFCXX4 REQUIRED netcdf-cxx4 IMPORTED_TARGET)
  link_libraries(PkgConfig::NETCDFCXX4)
endif (USE_NETCDF)

#strip leading whitespace from EXTRA_LIBS
string(REGEX REPLACE "^ " "" EXTRA_LIBS "${EXTRA_LIBS}")

set(EXECUTABLE_OUTPUT_PATH "bin")

message("FEASST_VERBOSE_LEVEL ${FEASST_VERBOSE_LEVEL}")
add_definitions("-DFEASST_VERBOSE_LEVEL_=${FEASST_VERBOSE_LEVEL}")
#unset(FEASST_VERBOSE_LEVEL CACHE)
add_definitions("-DFEASST_VERSION_=\"${FEASST_VERSION}\"")
add_definitions("-DFEASST_DIR_=\"${CMAKE_SOURCE_DIR}\"")

# OMP
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()
if (USE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    message("OpenMP found")
    include_directories(${OpenMP_C_INCLUDE_DIR})
    if(${OpenMP_CXX_INCLUDE_DIRS})
      include_directories( ${OpenMP_CXX_INCLUDE_DIRS} )
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endif()

if (EV)
  # prep feasst.h (and documentation)
  execute_process(
    COMMAND python3 ${CMAKE_SOURCE_DIR}/dev/tools/depend.py -s ${CMAKE_SOURCE_DIR} --update_doc 0 --check_headers ${USE_HEADER_CHECK}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    RESULT_VARIABLE EXIT_CODE
  )

  if (NOT ${EXIT_CODE} EQUAL 0)
      message(
              FATAL_ERROR
              "There was a problem generating feasst.h using the following command: \"python3 ../dev/tools/depend.py -s ../\" Ensure that python3 is available, and try to run this command manually to report any issues to the developers."
      )
  endif()
endif (EV)

# FEASST plugins
set (FEASST_PLUGIN_DIR "${CMAKE_SOURCE_DIR}/plugin")
include_directories (${FEASST_PLUGIN_DIR})
foreach (PLUGIN ${FEASST_PLUGINS})
  set(PLUGIN_DIR "${FEASST_PLUGIN_DIR}/${PLUGIN}")
  # create CMAKE_FEASST_PLUGIN_DIRS for Doxygen
  set(CMAKE_FEASST_PLUGIN_DIRS "${CMAKE_FEASST_PLUGIN_DIRS} ${PLUGIN_DIR}/include")
  file (GLOB src_cpp "${PLUGIN_DIR}/src/*.cpp")
  add_library (feasst${PLUGIN} SHARED ${src_cpp})
  target_link_libraries(feasst${PLUGIN} ${EXTRA_LIBS})
  if (OPENMP_FOUND)
    target_link_libraries(feasst${PLUGIN} OpenMP::OpenMP_CXX)
  endif()
endforeach (PLUGIN)

# HWH fix this: clang seems to require these dependencies but g++ does not... ?

# utils
list (FIND FEASST_PLUGINS "utils" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstutils feasstthreads)
endif()

# math
list (FIND FEASST_PLUGINS "math" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmath feasstutils)
endif()

# configuration
list (FIND FEASST_PLUGINS "configuration" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstconfiguration feasstmath)
endif()

# system
list (FIND FEASST_PLUGINS "system" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstsystem feasstconfiguration)
endif()

# opt_lj
list (FIND FEASST_PLUGINS "opt_lj" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstopt_lj feasstsystem)
endif()

# monte_carlo
list (FIND FEASST_PLUGINS "monte_carlo" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmonte_carlo feasstsystem)
endif()

# steppers
list (FIND FEASST_PLUGINS "steppers" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasststeppers feasstmonte_carlo)
endif()

# models
list (FIND FEASST_PLUGINS "models" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmodels feasstsystem feasstexample)
endif()

# mayer
list (FIND FEASST_PLUGINS "mayer" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmayer feasstmonte_carlo)
endif()

# patch
list (FIND FEASST_PLUGINS "patch" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstpatch feasstmonte_carlo feasstmayer)
endif()

# example
list (FIND FEASST_PLUGINS "example" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstexample feasstsystem feasstmonte_carlo feasststeppers)
endif()

# shape
list (FIND FEASST_PLUGINS "shape" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstshape feasstmath)
endif()

# confinement
list (FIND FEASST_PLUGINS "confinement" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstconfinement feasstmonte_carlo feasstshape)
endif()

# cluster
list (FIND FEASST_PLUGINS "cluster" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstcluster feasstmonte_carlo)
endif()

# chain
list (FIND FEASST_PLUGINS "chain" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstchain feasstmonte_carlo feasstcluster feasstgibbs)
endif()

# flat_histogram
list (FIND FEASST_PLUGINS "flat_histogram" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstflat_histogram feasstmonte_carlo feasststeppers feasstthreads)
endif()

# charge
list (FIND FEASST_PLUGINS "charge" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstcharge feasstmonte_carlo)
endif()

# morph
list (FIND FEASST_PLUGINS "morph" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmorph feasstchain feasstflat_histogram feasstcharge)
endif()

# prefetch
list (FIND FEASST_PLUGINS "prefetch" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstprefetch feasstmonte_carlo feasstthreads)
endif()

# egce
list (FIND FEASST_PLUGINS "egce" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstegce feasstflat_histogram)
endif()

# beta_expanded
list (FIND FEASST_PLUGINS "beta_expanded" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstbeta_expanded feasstflat_histogram)
endif()

# aniso
list (FIND FEASST_PLUGINS "aniso" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstaniso feasstmonte_carlo)
endif()

# gibbs
list (FIND FEASST_PLUGINS "gibbs" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstgibbs feasstmonte_carlo)
endif()

# model_expanded
list (FIND FEASST_PLUGINS "model_expanded" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmodel_expanded feasstflat_histogram)
endif()

# server
list (FIND FEASST_PLUGINS "server" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstserver feasstmonte_carlo)
endif()

# mpi
list (FIND FEASST_PLUGINS "mpi" _index)
if (${_index} GREATER -1)
  find_package(MPI REQUIRED)
  target_link_libraries(feasstmpi feasstmonte_carlo MPI::MPI_CXX)
endif()

# fftw
list (FIND FEASST_PLUGINS "fftw" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstfftw feasstmonte_carlo)
endif()

# netcdf
list (FIND FEASST_PLUGINS "netcdf" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstnetcdf feasstmonte_carlo)
endif()

# xtc
list (FIND FEASST_PLUGINS "xtc" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstxtc feasstsystem)
endif()

# CCACHE
if (USE_CCACHE)
  find_program(CCACHE_PROGRAM ccache)
  if(CCACHE_PROGRAM)
    message("CCACHE FOUND")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  endif(CCACHE_PROGRAM)
endif(USE_CCACHE)

# SPHINX
if (USE_SPHINX)
  find_package(Doxygen REQUIRED)
  if(DOXYGEN_FOUND)
    #message("Doxygen found")
    configure_file(${CMAKE_SOURCE_DIR}/dev/sphinx/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doxygen
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
  endif(DOXYGEN_FOUND)

  find_package(Sphinx REQUIRED)

  # set(SPHINX_THEME "nature")
  set(SPHINX_THEME "sphinx_rtd_theme")

  # configured documentation tools and intermediate build results
  set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")

  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  # LATEX output directory
  set(SPHINX_LATEX_DIR "${CMAKE_CURRENT_BINARY_DIR}/latex")

  configure_file(
    "${CMAKE_SOURCE_DIR}/dev/sphinx/conf.py.in"
    "${BINARY_BUILD_DIR}/conf.py"
    @ONLY)

  configure_file(
    "${CMAKE_SOURCE_DIR}/dev/sphinx/README.rst.in"
    "${CMAKE_SOURCE_DIR}/README.rst"
    @ONLY)

  configure_file(
    "${CMAKE_SOURCE_DIR}/dev/sphinx/README_html.rst.in"
    "${CMAKE_SOURCE_DIR}/README_html.rst"
    @ONLY)

  # copy favicon, leave notice and extended layout
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_static/favicon.ico ${BINARY_BUILD_DIR}/favicon.ico COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_static/leave_notice.js ${BINARY_BUILD_DIR}/_static/leave_notice.js COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_static/redirect_alert.js ${BINARY_BUILD_DIR}/_static/redirect_alert.js COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/layout.html ${BINARY_BUILD_DIR}/_templates/layout.html COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/footer.html ${BINARY_BUILD_DIR}/_templates/footer.html COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/breadholder.html ${BINARY_BUILD_DIR}/_templates/breadholder.html COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/breadcrumbs.html ${BINARY_BUILD_DIR}/_templates/breadcrumbs.html COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/boilerplate-header.html ${BINARY_BUILD_DIR}/_templates/boilerplate-header.html COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/boilerplate-footer.html ${BINARY_BUILD_DIR}/_templates/boilerplate-footer.html COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_templates/searchbox.html ${BINARY_BUILD_DIR}/_templates/searchbox.html COPYONLY)

  add_custom_target(html ALL
    ${SPHINX_EXECUTABLE}
        -q -b html
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "../"
        "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx")

  add_dependencies(html doxygen)

  add_custom_target(latex ALL
    ${SPHINX_EXECUTABLE}
        -q -b latex
        -c "${BINARY_BUILD_DIR}"
        -d "${SPHINX_CACHE_DIR}"
        "../"
        "${SPHINX_LATEX_DIR}"
    COMMENT "Building LATEX documentation with Sphinx")

  add_dependencies(latex doxygen)

  add_custom_target(pdf ALL
    ${CMAKE_COMMAND} -E chdir "${SPHINX_LATEX_DIR}" make
    COMMENT "Building PDF from LATEX documentation with Sphinx")

  add_dependencies(pdf latex)

endif (USE_SPHINX)

# GTEST
if (USE_GTEST)
  message("USING GTEST")
  enable_testing()
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/tmp")
  file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/build-gtest")
  # Download and unpack googletest at configure time
  configure_file(${CMAKE_SOURCE_DIR}/dev/cmake/gtest.cmake ${CMAKE_SOURCE_DIR}/build-gtest/googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build-gtest/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/build-gtest/googletest-download )
  if(result)
    message("Build step for googletest failed: ${result}")
  endif()

  # Prevent overriding the parent project's compiler/linker
  # settings on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_SOURCE_DIR}/build-gtest/googletest-src
                   ${CMAKE_SOURCE_DIR}/build-gtest/googletest-build
                   EXCLUDE_FROM_ALL)

  # The gtest/gtest_main targets carry header search path
  # dependencies automatically when using CMake 2.8.11 or
  # later. Otherwise we have to add them here ourselves.
  if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
  endif()

  # Add test source files from each lib
  foreach (PLUGIN ${FEASST_PLUGINS})
    file(GLOB lib_test_src "${FEASST_PLUGIN_DIR}/${PLUGIN}/test/*.cpp")
    set(test_src ${test_src} ${lib_test_src})
  endforeach (PLUGIN)

  # Now simply link against gtest or gtest_main as needed. Eg
  add_executable(unittest ${test_src})
  target_link_libraries(unittest gtest_main ${EXTRA_LIBS})

  foreach (PLUGIN ${FEASST_PLUGINS})
    target_link_libraries(unittest feasst${PLUGIN})
  endforeach (PLUGIN)
  add_test(NAME test COMMAND ./bin/unittest --gtest_shuffle --gtest_filter=-*LONG)

endif(USE_GTEST)

# build feasst library
add_library (feasstlib INTERFACE)
foreach (PLUGIN ${FEASST_PLUGINS})
  target_link_libraries(feasstlib INTERFACE feasst${PLUGIN})
endforeach (PLUGIN)
target_link_libraries(feasstlib INTERFACE ${EXTRA_LIBS})
#list (FIND FEASST_PLUGINS "mpi" _index)
#if (${_index} GREATER -1)
#  message("linking mpi")
#  target_link_libraries(feasstlib INTERFACE MPI::MPI_CXX)
#endif()

# install library inside the binary "build" directory
#file (WRITE ${CMAKE_BINARY_DIR}/feasst.h "")
foreach (PLUGIN ${FEASST_PLUGINS})
  target_include_directories(feasst${PLUGIN} INTERFACE
      $<BUILD_INTERFACE:${FEASST_PLUGIN_DIR}/${PLUGIN}/include>
      $<INSTALL_INTERFACE:include/${PLUGIN}/include>)
  file(GLOB feasst_headers "${FEASST_PLUGIN_DIR}/${PLUGIN}/include/*.h")
  install(FILES ${feasst_headers} DESTINATION include/${PLUGIN}/include)
  file(GLOB feasst_headers RELATIVE "${FEASST_PLUGIN_DIR}" "${FEASST_PLUGIN_DIR}/${PLUGIN}/include/*.h")
  # generate header feasst.h
#  foreach (HEADER ${feasst_headers})
#    file (APPEND ${CMAKE_BINARY_DIR}/feasst.h "#include \"${HEADER}\"\n")
#  endforeach (HEADER)
endforeach (PLUGIN)
#message(${CMAKE_SOURCE_DIR}/plugin/feasst/include/feasst.h)
#message(${CMAKE_CURRENT_BINARY_DIR}/feasst.h)
configure_file(
  ${CMAKE_SOURCE_DIR}/plugin/feasst/include/feasst.h ${CMAKE_CURRENT_BINARY_DIR}/feasst.h COPYONLY)
foreach (PLUGIN ${FEASST_PLUGINS})
  message(feasst${PLUGIN})
  install(TARGETS feasst${PLUGIN}
      EXPORT feasst-targets
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib
      RUNTIME DESTINATION bin)
endforeach (PLUGIN)
install(TARGETS feasstlib
    EXPORT feasst-targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
install(EXPORT feasst-targets
    FILE feasst-config.cmake
    DESTINATION lib/cmake/feasst)
install(FILES ${CMAKE_BINARY_DIR}/feasst.h
    DESTINATION include)

# make feasst text interpreter
list (FIND FEASST_PLUGINS "monte_carlo" _index)
if (${_index} GREATER -1)
  add_executable (fst ${CMAKE_SOURCE_DIR}/plugin/feasst/src/feasst.cpp)
  target_link_libraries (fst LINK_PUBLIC feasstlib)
  add_executable (rst ${CMAKE_SOURCE_DIR}/plugin/feasst/src/restart.cpp)
  target_link_libraries (rst LINK_PUBLIC feasstlib)
endif()

# PYBIND11
if (USE_PYBIND11)
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  find_package(pybind11 CONFIG REQUIRED)
  python3_add_library(_core MODULE ${CMAKE_SOURCE_DIR}/python/feasst.cpp WITH_SOABI)
  target_link_libraries(_core PRIVATE pybind11::headers)
  target_link_libraries (_core LINK_PUBLIC feasstlib)
  target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
  install(TARGETS _core DESTINATION feasst)
endif (USE_PYBIND11)
