#!/usr/bin/env bash
if [ -z $1 ]; then
  echo "Usage: ./run_py.sh [name of python script].py [optional arguments]"
  exit
fi

# Obtain the current directory of the script.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

PYTHONPATH=$PYTHONPATH:${DIR}/../build python3 $@
