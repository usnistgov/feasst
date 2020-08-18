#!/bin/bash

mkdir -p build
cd build
python3 -m venv feasst_test_env
source feasst_test_env/bin/activate
pip install jupyter matplotlib pandas scipy
cmake -DUSE_GTEST=ON -DUSE_SWIG=ON ..
make _feasst -j8
make install -j8
echo "" > summary_long.log
echo "" > summary.log

valgrind ./bin/unittest --gtest_filter=-*LONG* >> summary_long.log 2>&1
echo "********** valgrind **********" >> summary.log
tail -1 summary_long.log >> summary.log

python ../py/test.py >> summary_long.log 2>&1
echo "********** pytest **********" >> summary.log
tail -1 summary_long.log >> summary.log

echo "********** pytutorial **********" >> summary.log
python ../dev/tools/run_tutorials.py >> summary_long.log 2>&1
tail -1 tutorial_failures.txt >> summary.log

