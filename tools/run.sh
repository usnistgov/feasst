#!/usr/bin/env bash
# Obtain the current directory of the script.
# It should be in tools, where the compilation script is also located.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# check that the file extension is py, c, cc or cpp
filename=$(basename "$1")
extension="${filename##*.}"
##filename="${filename%.*}"
if [[ $extension == "py" ]]; then
  $DIR/run_py.sh $@
else
  if [[ $extension != "cc" ]] && [[ $extension != "c" ]] && [[ $extension != "cpp" ]]; then
    echo "Unrecognized file extension: \"$extension\". Expecting py, c, cc, or cpp"
    exit
  else
    $DIR/run_cc.sh $@
  fi
fi
