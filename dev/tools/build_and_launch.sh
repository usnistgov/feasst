#!/bin/bash
mkdir -p build
cd build
cmake ..
make install -j12

python3 -m venv feasst_test_env
source feasst_test_env/bin/activate
pip install pyfeasst numpy jupyter matplotlib pandas scipy

echo "" > summary.log
echo "" > summary_long.log

#tail -1 tutorial_failures.txt >> summary.log
echo "********** launch py **********" >> summary.log
python ../dev/tools/lnch_tutorials.py >> summary_long.log 2>&1
grep Error summary_long.log >> summary.log
cat launch_failures.txt >> summary.log

