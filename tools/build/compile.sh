file=$(basename $0 .sh)
cmake -DUSE_GTEST=ON -DUSE_GTEST_TUTORIALS=ON -DUSE_CCACHE=ON -DUSE_SWIG=ON .. > $file.log 2>&1
make -j 8 >> $file.log 2>&1
make -j 8 unittest >> $file.log 2>&1
echo "**** $file ****" >> summary.log
tail -1 $file.log >> summary.log
echo "" >> summary.log
