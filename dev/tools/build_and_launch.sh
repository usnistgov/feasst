#!/bin/bash
mkdir -p build
cd build
cmake ..
make install -j12
echo "" > summary.log
echo "" > summary_long.log

#tail -1 tutorial_failures.txt >> summary.log
echo "********** launch py **********" >> summary.log
python ../dev/tools/lnch_tutorials.py >> summary_long.log 2>&1
cat launch_failures.txt >> summary.log

