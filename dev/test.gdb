## Usage: gdb -batch -x filename executable
catch throw
#r < script2.txt
#r < pore000_fstin.txt
r --gtest_filter=MonteCarlo.gibbs_ensemble
bt
