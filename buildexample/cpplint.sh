file=$(basename $0 .sh)
/usr/bin/python ../tools/cpplint.py ../src/* > $file.log 2>&1
echo "**** $file ****" >> summary.log
tail -4 $file.log >> summary.log
echo "" >> summary.log
