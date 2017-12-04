#!/usr/bin/env bash
DES=testcase
mkdir -p $DES

# generate aggregated c++ testcases
for file in `find ../testcase -name test.cc`; do
  newfile=$DES/`echo $file | sed 's/\.\.\/testcase\///' | sed 's/\//_/g'`
  ../tools/cc2gtest.sh $file $newfile
done

# generate aggregated python testcases
pyfile=$DES/test.py
echo "" > $pyfile
for file in `find ../testcase -name test.py`; do
  head -n -3 $file >> $pyfile
done
cat <<-EOF >> $pyfile
if __name__ == "__main__":
    unittest.main()
EOF

# generate / filter RST files for documentation
for file in `find ../testcase -name *README.rst | grep -v "bak"`; do
  DIR=$(dirname "${file}" | sed 's/\.\.\///' | sed 's/\//\\\//g')
  newfile=$DES/`echo $file | sed 's/\.\.\/testcase\///' | sed 's/\//_/g'`
  cat $file | sed "s/AUTO_GEN_DIR/$DIR/" > $newfile
done

# finally, test that there are no broken HTML links
for file in `ls ../src/*cc ../src/*h`; do
  rstfile=`grep "<a href=" $file | sed 's/html/rst/' | sed 's/\(.*\)<a href=\"\(.*\)\">\(.*\)<\/a>\(.*\)/\2/'`
  if [ ! -f $rstfile ]; then
    echo "ERROR: Broken HTML link in file $file which references $rstfile"
    exit
  fi
done

###################
# begin api.rst gen
###################

# Note that the following classes must be defined in a very specific order
baseClasses="accumulator analyze barrier base criteria custom_exception histogram mc pair random shape space table trial"
nonClasses="functions ui_abbreviated physical_constants"

#####################################
echo "Generating documentation *.rst"
#####################################

printRSTHeader()
{
  className=`grep "class " $header | grep "{" | head -1 | awk '{print $2}'`
  printfile=${className}.rst
  echo "   $className" >> $tocfile
  echo "$className" > $printfile
  echo "=====================================================" >> $printfile
  echo "" >> $printfile
  if [ -z "$nonClassFlag" ]; then
    echo ".. doxygenclass:: $className" >> $printfile
    echo "   :project: FEASST" >> $printfile
    echo "   :members:" >> $printfile
  else
    echo ".. doxygenfile:: `basename $header`" >> $printfile
    echo "   :project: FEASST" >> $printfile
  fi
}

printDoc()
{
  arr=(`ls -vd ../src/${class}_*.h 2> /dev/null`)
  # if array finds derived classes, set up toc
  header="../src/${class}.h"
  className=`grep "class " $header | grep "{" | head -1 | awk '{print $2}'`
  tocfile=${class}toc.rst
  echo "$className Classes" > $tocfile
  echo "=======================================================" >> $tocfile
  echo "" >> $tocfile
  echo ".. toctree::" >> $tocfile
  printRSTHeader
  if (( ${#arr[@]} > 0 )); then
    for header in ${arr[@]}; do
      printRSTHeader
    done
    echo "   ${class}toc" >> $apifile
  else
    echo "   $className" >> $apifile
  fi
}

# begin printing the api file, which needs to know whether to include one
# class, or a toc for may derived classes
apifile=api.rst
cat <<-EOF > $apifile
***********************************
Application program interface (API)
***********************************

This section describes the API which is available for use via C++ or python.

.. toctree::
EOF

# print the rst for each class and also add to the api
nonClassFlag=""
for class in $baseClasses; do
  printDoc
done
nonClassFlag="True"
for class in $nonClasses; do
  printDoc
done

