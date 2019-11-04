mkdir -p build
cd build
cmake -DUSE_GTEST=ON -DUSE_SWIG=ON ..
make -j 4
echo "" > summary_long.log
echo "" > summary.log

valgrind ./bin/unittest >> summary_long.log 2>&1
echo "********** valgrind **********" >> summary.log
tail -1 summary_long.log >> summary.log

python ../py/test.py >> summary_long.log 2>&1
echo "********** pytest **********" >> summary.log
tail -1 summary_long.log >> summary.log

echo "********** pytutorial **********" >> summary.log
python ../dev/tools/run_tutorials.py >> summary_long.log 2>&1
tail -1 tutorial_failures.txt >> summary.log

