## Usage: gdb -batch -x filename executable
catch throw
#r < run.txt
r < trappe000_fstin.txt
#r < sal000_fstin.txt
#r --gtest_filter=MonteCarlo.TrialMorphExpanded
bt
