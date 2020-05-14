#!/usr/bin/env bash
if [ -z $1 ]; then
  echo "Usage: ./run_py.sh [name of python script].py [optional arguments]"
  exit
fi

# Obtain the current directory of the script.
# It should be in tools, where the compilation script is also located.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p tmp  # make a tmp directory for restart files
PYTHONPATH=$PYTHONPATH:${DIR}/../build python $@
