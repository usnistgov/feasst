#!/usr/bin/env bash
# This script is used to compile FEASST input files using the C++ interface
# Usage: /path/to/feasst/tools/compile.sh [program name] [build dirname: opt]

prog=$1
builddir=$2
if [ -z $builddir ]; then
  builddir="build"
fi

# silence pushd
pushd() {
    command pushd "$@" > /dev/null
}

# silence popd
popd() {
    command popd "$@" > /dev/null
}

# Obtain the current directory of the script.
# It should be in tools
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cp ${prog}.cc $DIR/../drivers/main.cc
pushd $DIR/../$builddir
  echo "Using $0 to build $1.cc in "`pwd`
  make -j main || exit 1
popd
cp $DIR/../$builddir/bin/main $prog

