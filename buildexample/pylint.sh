file=$(basename $0 .sh)
pylint ../src/*.py > $file.log 2>&1
pylint ../src/swigtest/*.py > $file.log 2>&1
echo "**** $file ****" >> summary.log
tail -4 $file.log >> summary.log
echo "" >> summary.log
