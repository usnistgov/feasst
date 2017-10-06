file=$(basename $0 .sh)
make coverage > $file.log 2>&1
echo "**** $file ****" >> summary.log
tail -4 $file.log >> summary.log
echo "" >> summary.log
