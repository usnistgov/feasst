file=$(basename $0 .sh)
../tools/run.sh tutorial/test.py > $file.log 2>&1
echo "**** $file ****" >> summary.log
grep "OK" $file.log >> summary.log
grep "FAIL" $file.log >> summary.log
echo "" >> summary.log
