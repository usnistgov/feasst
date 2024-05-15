#!/bin/bash
export OMP_NUM_THREADS=4
feasst_dir=`pwd`
mkdir -p build
cd build
python3 -m venv feasst_test_env
source feasst_test_env/bin/activate
python3 -m pip install --upgrade pip
pip install ../pyfeasst numpy jupyter matplotlib pandas scipy wheel biopandas
module load mpi/openmpi-x86_64 # sudo dnf install openmpi-devel
cmake -DUSE_GTEST=ON -DUSE_MPI=ON .. #-DUSE_SWIG=ON
#make feasst -j24
make install -j24
echo "" > summary_long.log
echo "" > summary.log

valgrind ./bin/unittest --gtest_color=no --gtest_filter=-*LONG* >> summary_long.log 2>&1
echo "********** valgrind **********" >> summary.log
grep FAIL summary_long.log >> summary.log
grep Throw summary_long.log >> summary.log
tail -1 summary_long.log >> summary.log

~/software/gtest-parallel/gtest-parallel bin/unittest --gtest_color=no -d tmp >> summary_long_allt.log 2>&1
cat summary_long_allt.log >> summary_long.log
echo "********** all gtest **********" >> summary.log
grep FAIL summary_long_allt.log >> summary.log
grep Throw summary_long_allt.log >> summary.log

pushd ../pyfeasst/src/pyfeasst
  for f in *.py; do
    python "$f" >> summary_pyfeasst.log
  done
popd
echo "********** pyfeasst **********" >> summary.log
tail -1 ../pyfeasst/src/pyfeasst/summary_pyfeasst.log >> summary.log

#python ../py/test.py >> summary_long.log 2>&1
#echo "********** pytest **********" >> summary.log
#tail -1 summary_long.log >> summary.log

echo "********** pytutorial **********" >> summary.log
python ../dev/tools/run_tutorials.py --feasst_install $feasst_dir >> summary_long.log 2>&1
for fl in `find ../ -name 'tutorial_failures.txt'`; do
  cat $fl >> summary.log
done

#tail -1 tutorial_failures.txt >> summary.log
#echo "********** launch py **********" >> summary.log
#python ../dev/tools/lnch_tutorials.py >> summary_long_long.log 2>&1
#tail -1 launch_failures.txt >> summary.log

