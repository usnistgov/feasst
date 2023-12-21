cmake_minimum_required (VERSION 3.5)
project (FEASST)
set (CMAKE_CXX_STANDARD 11)
cmake_policy(SET CMP0042 NEW) # removes warning about macosx_rpath
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -g -pthread")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -Wall -pedantic -g -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -g -pthread")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pg -g") #profile

# manually list the plugins
# Note: for developement and testing, reducing plugins decreases compile time
#set(FEASST_PLUGINS "threads")
#set(FEASST_PLUGINS "threads;utils")
#set(FEASST_PLUGINS "threads;utils;math")
#set(FEASST_PLUGINS "threads;utils;math;configuration")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;models;monte_carlo")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;aniso")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;cluster")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;prefetch;charge;flat_histogram;cluster")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;prefetch;charge;flat_histogram;cluster;beta_expanded")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;prefetch;charge;flat_histogram;cluster;beta_expanded;chain")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;prefetch;charge;flat_histogram;cluster;beta_expanded;chain;morph")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;mayer;patch;example;shape;confinement;flat_histogram;charge;cluster;prefetch;egce;aniso")
set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;mayer;patch;example;shape;confinement;chain;flat_histogram;charge;morph;cluster;prefetch;egce;beta_expanded;aniso;gibbs")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;mayer;patch;example;shape;confinement;chain;flat_histogram;charge;morph;cluster;prefetch;egce;beta_expanded;aniso;gibbs;fftw")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;mayer;patch;example;shape;confinement;chain;flat_histogram;charge;morph;cluster;prefetch;egce;beta_expanded;aniso;netcdf")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;mayer;patch;example;shape;confinement;chain;flat_histogram;charge;morph;cluster;prefetch;egce;beta_expanded;aniso;fftw")
#set(FEASST_PLUGINS "threads;utils;math;configuration;system;opt_lj;steppers;monte_carlo;models;patch;mayer;example;shape;confinement;chain;flat_histogram;charge;morph;cluster;prefetch;egce;beta_expanded;aniso;xtc")
# Note, if you change the plugin list and want to use the python interface,
# you also need to run "../py/run.sh ../dev/tools/depend.py -s ../"

# set path to XDR if using xtc plugin
set(XDRFILE_DIR "$ENV{HOME}/software/xdrfile-1.1b/build")

option(USE_GTEST "Use gtest" OFF)
option(USE_CCACHE "Use ccache to speed up builds" ON)
option(USE_OPENMP "Use OpenMP" ON)
option(USE_SPHINX "Use SPHINX for documentation" OFF)
option(USE_SWIG "Use SWIG for python interface" OFF)
option(USE_GCOV "Use Coverage" OFF)
option(USE_FFTW "Use FFTW" OFF)
set(FFTW_DIR "$ENV{HOME}/software/fftw-3.3.10/build/")
option(USE_NETCDF "Use NetCDF" OFF)
set(NETCDF_DIR "$ENV{HOME}/local")

# If you would like to use the python interface, you may need to specify the
# python library and include path variables below. If not specified then
# CMake attempts to find the libraries on its own.
option(SET_PYTHON_PATH "Manually define PYTHON_INCLUDE_DIR and PYTHON_LIBRARY" OFF)

# Users typically would not edit below this comment

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
execute_process(
  COMMAND git describe --abbrev=0
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH_SHORT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)


add_definitions("-DFEASST_VERSION_=\"${GIT_COMMIT_HASH}-${GIT_BRANCH}\"")
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
  target_link_libraries(feasstmodels feasstsystem)
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

  # FEASST version
  set(FEASST_VERSION ${GIT_COMMIT_HASH})

  configure_file(
    "${CMAKE_SOURCE_DIR}/dev/sphinx/conf.py.in"
    "${BINARY_BUILD_DIR}/conf.py"
    @ONLY)

  # copy favicon, leave notice and extended layout
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_static/favicon.ico ${BINARY_BUILD_DIR}/favicon.ico COPYONLY)
  configure_file(
    ${CMAKE_SOURCE_DIR}/dev/sphinx/_build/_static/leave_notice.js ${BINARY_BUILD_DIR}/_static/leave_notice.js COPYONLY)
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

# prep feasst.(i/h) (as well as documentation)
execute_process(
  COMMAND python3 ../py/depend.py -s ../ --update_doc 0
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  RESULT_VARIABLE EXIT_CODE
)

if (NOT ${EXIT_CODE} EQUAL 0)
    message(
            FATAL_ERROR
            "There was a problem generating feasst.h using the following command: \"python3 ../py/depend.py -s ../\" Ensure that python3 is available, and try to run this command manually to report any issues to the developers."
    )
endif()

# build feasst library
add_library (feasstlib INTERFACE)
foreach (PLUGIN ${FEASST_PLUGINS})
  target_link_libraries(feasstlib INTERFACE feasst${PLUGIN})
endforeach (PLUGIN)
target_link_libraries(feasstlib INTERFACE ${EXTRA_LIBS})

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

# SWIG
if (USE_SWIG)
  if(POLICY CMP0078)
    cmake_policy(SET CMP0078 NEW)
  endif()
  if(POLICY CMP0086)
    cmake_policy(SET CMP0086 NEW)
  endif()
  if (SET_PYTHON_PATH)
  else (SET_PYTHON_PATH)
    FIND_PACKAGE(PythonLibs 3)
    if (PythonLibs_FOUND)
    else (PythonLibs_FOUND)
      # Try default anaconda home install.
      set(PYTHON_INCLUDE_DIR "$ENV{HOME}/anaconda3/include/python3.7m/")
      set(PYTHON_LIBRARY "$ENV{HOME}/anaconda3/lib/libpython3.7m.so")
    endif (PythonLibs_FOUND)
  endif (SET_PYTHON_PATH)
  message("PYTHON_INCLUDE_DIR: ${PYTHON_INCLUDE_DIR}")
  message("PYTHON_LIBRARY: ${PYTHON_LIBRARY}")

  FIND_PACKAGE(SWIG)
  if (SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_DIR})
    #INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    get_filename_component(barename "${CMAKE_CURRENT_BINARY_DIR}" NAME)
    add_custom_target(build-time-make-directory ALL
      COMMAND ${CMAKE_COMMAND} -E make_directory ${barename})
    SET(CMAKE_SWIG_FLAGS "")

    SET_SOURCE_FILES_PROPERTIES("${CMAKE_SOURCE_DIR}/py/feasst.i" PROPERTIES CPLUSPLUS ON)
    #SET_SOURCE_FILES_PROPERTIES("${CMAKE_CURRENT_BINARY_DIR}/feasst.i" PROPERTIES SWIG_FLAGS "-includeall")
    if (${CMAKE_VERSION} VERSION_LESS "3.10.0")
      SWIG_ADD_MODULE(feasst python "${CMAKE_SOURCE_DIR}/py/feasst.i")
    else ()
      SWIG_ADD_LIBRARY(feasst LANGUAGE python SOURCES "${CMAKE_SOURCE_DIR}/py/feasst.i")
    endif ()
    SWIG_LINK_LIBRARIES(feasst ${EXTRA_LIBS})
    SWIG_LINK_LIBRARIES(feasst feasstlib)
    if (APPLE)
      if (${CMAKE_VERSION} VERSION_LESS "3.15.0")
        set_target_properties(_feasst PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
      else ()
        set_target_properties(feasst PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
      endif ()
    endif (APPLE)
    if (${CMAKE_VERSION} VERSION_LESS "3.15.0")
      add_dependencies(_feasst build-time-make-directory)
    else ()
      add_dependencies(feasst build-time-make-directory)
    endif ()

    # install using setup.py
    set(OUTPUT    "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp")
    set(SETUP_PY  "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    configure_file("${CMAKE_SOURCE_DIR}/dev/tools/setup.py.in" ${SETUP_PY})
#    add_custom_command(OUTPUT ${OUTPUT}
#                       COMMAND python3 setup.py build
#                       COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT})
#    add_custom_target(target ALL DEPENDS ${OUTPUT})
    install(CODE "execute_process(COMMAND pip install .)")
    #install(CODE "execute_process(COMMAND python3 ${SETUP_PY} install)")

    add_test(NAME pytest COMMAND python3 ${CMAKE_SOURCE_DIR}/py/test.py)
    add_test(NAME pytutorial COMMAND python3 ../dev/tools/run_tutorials.py)
  else (SWIG_FOUND)
    message("SWIG NOT FOUND")
  endif (SWIG_FOUND)
endif(USE_SWIG)
