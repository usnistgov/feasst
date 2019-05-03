cmake_minimum_required (VERSION 2.6)
project (FEASST)
set (CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -g -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -g -pthread")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pg -g") #profile

# manually list the plugins
#set(FEASST_PLUGINS "core;")
#set(FEASST_PLUGINS "core;ewald;example;")
#set(FEASST_PLUGINS "core;ewald;example;confinement;xtc;")
set(FEASST_PLUGINS "core;ewald;example;confinement;chain;models;patch;flat_histogram;mayer;steppers")
message(${FEASST_PLUGINS})
# set path to XDR if using xtc plugin
set(XDRFILE_DIR "$ENV{HOME}/software/xdrfile-1.1b/build")

# set the install location (default: in the binary build directory)"
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

option(USE_GTEST "Use gtest" ON)
option(USE_CCACHE "Use ccache to speed up builds" ON)
option(USE_SPHINX "Use SPHINX for documentation" OFF)
option(USE_SWIG "Use SWIG for python interface" OFF)
option(USE_GCOV "Use Coverage" OFF)

# If you would like to use the python interface, you may need to specify the
# python library and include path variables below. If not specified then
# CMake attempts to find the libraries on its own.
option(SET_PYTHON_PATH "Manually define PYTHON_INCLUDE_PATH and PYTHON_LIBRARIES" OFF)

# Users typically would not edit below this comment

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
  if(CMAKE_COMPILER_IS_GNUCXX)
    include(CodeCoverage)
    setup_target_for_coverage(coverage unittest coverage)
  endif()
endif (USE_GCOV)

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

add_definitions("-DFEASST_VERSION_=\"${GIT_COMMIT_HASH} ${GIT_BRANCH}\"")
add_definitions("-DFEASST_DIR_=\"${CMAKE_SOURCE_DIR}\"")

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
endforeach (PLUGIN)

# HWH fix this: clang seems to require these dependencies but g++ does not... ?
# XTC
list (FIND FEASST_PLUGINS "xtc" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstxtc feasstcore)
endif()

# Ewald
list (FIND FEASST_PLUGINS "ewald" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstewald feasstcore)
endif()

# Example
list (FIND FEASST_PLUGINS "example" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstexample feasstcore)
endif()

# Confinement
list (FIND FEASST_PLUGINS "confinement" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstconfinement feasstcore)
endif()

# patch
list (FIND FEASST_PLUGINS "patch" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstpatch feasstcore)
endif()

# models
list (FIND FEASST_PLUGINS "models" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmodels feasstcore)
endif()

# chain
list (FIND FEASST_PLUGINS "chain" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstchain feasstcore)
endif()

# flat_histogram
list (FIND FEASST_PLUGINS "flat_histogram" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstflat_histogram feasstcore)
endif()

# mayer
list (FIND FEASST_PLUGINS "mayer" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasstmayer feasstcore)
endif()

# steppers
list (FIND FEASST_PLUGINS "steppers" _index)
if (${_index} GREATER -1)
  target_link_libraries(feasststeppers feasstcore)
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
                   ${CMAKE_SOURCE_DIR}/build-gtest/googletest-build)

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
  add_test(NAME test COMMAND ./bin/unittest --gtest_shuffle)

endif(USE_GTEST)

# build feasst library
add_library (feasst INTERFACE)
foreach (PLUGIN ${FEASST_PLUGINS})
  target_link_libraries(feasst INTERFACE feasst${PLUGIN})
endforeach (PLUGIN)
target_link_libraries(feasst INTERFACE ${EXTRA_LIBS})

# install library inside the binary "build" directory
file (WRITE ${CMAKE_BINARY_DIR}/feasst.h "")
foreach (PLUGIN ${FEASST_PLUGINS})
  target_include_directories(feasst${PLUGIN} INTERFACE
      $<BUILD_INTERFACE:${FEASST_PLUGIN_DIR}/${PLUGIN}/include>
      $<INSTALL_INTERFACE:include/${PLUGIN}/include>)
  file(GLOB feasst_headers "${FEASST_PLUGIN_DIR}/${PLUGIN}/include/*.h")
  install(FILES ${feasst_headers} DESTINATION include/${PLUGIN}/include)
  file(GLOB feasst_headers RELATIVE "${FEASST_PLUGIN_DIR}" "${FEASST_PLUGIN_DIR}/${PLUGIN}/include/*.h")
  # generate header feasst.h
  foreach (HEADER ${feasst_headers})
    file (APPEND ${CMAKE_BINARY_DIR}/feasst.h "#include \"${HEADER}\"\n")
  endforeach (HEADER)
endforeach (PLUGIN)
foreach (PLUGIN ${FEASST_PLUGINS})
  message(feasst${PLUGIN})
  install(TARGETS feasst${PLUGIN}
      EXPORT feasst-targets
      ARCHIVE DESTINATION lib
      LIBRARY DESTINATION lib
      RUNTIME DESTINATION bin)
endforeach (PLUGIN)
install(TARGETS feasst
    EXPORT feasst-targets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
install(EXPORT feasst-targets
    FILE feasst-config.cmake
    DESTINATION lib/cmake/feasst)
install(FILES ${CMAKE_BINARY_DIR}/feasst.h
    DESTINATION include)

# SWIG
if (USE_SWIG)
  if (SET_PYTHON_PATH)
  else (SET_PYTHON_PATH)
    FIND_PACKAGE(PythonLibs 3)
    if (PythonLibs_FOUND)
    else (PythonLibs_FOUND)
      # Try default anaconda home install.
      set(PYTHON_INCLUDE_PATH "$ENV{HOME}/anaconda3/include/python3.7m/")
      set(PYTHON_LIBRARIES "$ENV{HOME}/anaconda3/lib/libpython3.7m.so")
    endif (PythonLibs_FOUND)
  endif (SET_PYTHON_PATH)
  message("PYTHON_INCLUDE_PATH: ${PYTHON_INCLUDE_PATH}")
  message("PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}")

  FIND_PACKAGE(SWIG)
  if (SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    #INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    get_filename_component(barename "${CMAKE_CURRENT_BINARY_DIR}" NAME)
    add_custom_target(build-time-make-directory ALL
      COMMAND ${CMAKE_COMMAND} -E make_directory ${barename})
    SET(CMAKE_SWIG_FLAGS "")

    SET_SOURCE_FILES_PROPERTIES("${CMAKE_SOURCE_DIR}/py/feasst.i" PROPERTIES CPLUSPLUS ON)
    #SET_SOURCE_FILES_PROPERTIES("${CMAKE_CURRENT_BINARY_DIR}/feasst.i" PROPERTIES SWIG_FLAGS "-includeall")
    SWIG_ADD_MODULE(feasst python "${CMAKE_SOURCE_DIR}/py/feasst.i")
    SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS})
    SWIG_LINK_LIBRARIES(feasst feasst)
    add_dependencies(_feasst build-time-make-directory)

    # put pyfeasst into build
    configure_file(${CMAKE_SOURCE_DIR}/py/pyfeasst.py ${CMAKE_SOURCE_DIR}/build/pyfeasst.py @ONLY)

    add_test(NAME pytest COMMAND ${CMAKE_SOURCE_DIR}/py/run.sh ${CMAKE_SOURCE_DIR}/py/test.py)
  else (SWIG_FOUND)
    message("SWIG NOT FOUND")
  endif (SWIG_FOUND)
endif(USE_SWIG)
