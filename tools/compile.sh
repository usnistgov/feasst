#!/bin/bash
# This script is used to compile FEASST input files using the C++ interface
# Usage: /path/to/feasst/tools/compile.sh [program name] [build dirname: opt]
# First, the script decides whether or not cmake was used to build FEASST
# If it finds no evidence of cmake, then it tries to use the Makefile in src
# The $prog.cc file is placed in the source and compiled.
# Then the 'main' executable is placed in the current directory

prog=$1
builddir=$2
if [ -z $builddir ]; then
  builddir="build"
fi

# if the build directory exists with a Makefile, use cmake to compile
if [ -d $FEASST_INSTALL_DIR_/$builddir ] &&
   [ -s $FEASST_INSTALL_DIR_/$builddir/Makefile ]; then
  cp ${prog}.cc $FEASST_INSTALL_DIR_/drivers/main.cc
  pushd $FEASST_INSTALL_DIR_/$builddir
    echo "$0 is building in "`pwd`
    make main
  popd
  cp $FEASST_INSTALL_DIR_/$builddir/bin/main $prog

# otherwise, attempt to use the Makefile in src
else
  source=$FEASST_INSTALL_DIR_/src
  cp ${prog}.cc $source/main.cc
  pushd $source
    echo "$0 is building in "`pwd`
    make main
  popd
  cp $FEASST_INSTALL_DIR_/src/main $prog
fi
