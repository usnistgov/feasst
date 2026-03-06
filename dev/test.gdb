## Usage: gdb -batch -x filename executable
catch throw
#r < script2.txt
#r < test
#r < tr6d000_fstin.txt
r --gtest_filter=MonteCarlo.SineSlabTable_LONG
bt
