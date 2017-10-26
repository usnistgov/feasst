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

