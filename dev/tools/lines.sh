# Obtain the current directory of the script.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# all files
find $DIR/../../plugin/ -name '*.cpp' -o -name '*.h' | xargs wc -l | sort -n | tail -1 | awk '{print $1, "total"}'

# header files
find $DIR/../../plugin/ -name '*.cpp' -path '*/test/*' -o -name '*.h' -path '*/test/*' | xargs wc -l | sort -n | tail -1 | awk '{print $1, "tests"}'


