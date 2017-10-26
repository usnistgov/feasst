#!/usr/bin/env bash
if [ -z $1 ]; then
  echo "Usage: ./run.sh [name of C++ script with \"int main(\"].cc [optional arguments]"
  exit
fi

# check that the file extension is ".c", ".cc" or ".cpp"
filename=$(basename "$1")
extension="${filename##*.}"
##filename="${filename%.*}"
if [[ $extension != "cc" ]] && [[ $extension != "c" ]] && [[ $extension != "cpp" ]]; then
  echo "Unrecognized file extension: \"$extension\". Expecting c, cc, or cpp"
  exit
fi

# remove the file extension
prog=`echo $1 | sed 's/\.'"$extension"'//'`

# Obtain the current directory of the script.
# It should be in tools, where the compilation script is also located.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Warning: The next two lines MUST be run in this order to catch compile failure
$DIR/compile.sh $prog
FAILCODE=$?

if [ $FAILCODE -eq 0 ]; then
  mkdir -p tmp #directory for checkpoint files
  # remove first argument and send this to the program
  args=`echo "$@ " | sed 's/^[^ ]* //'`
  ./$prog $args
else
  echo "ERROR in $0: Compilation failed"
fi
