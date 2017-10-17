DES=testcase
mkdir -p $DES
for file in `find ../testcase/ -name *cc | grep -v "bak"`; do
  newfile=$DES/`echo $file | sed 's/\.\.\/testcase\///' | sed 's/\//_/g'`
  ../tools/cc2gtest.sh $file $newfile
  #echo $file
  #echo $newfile
done

