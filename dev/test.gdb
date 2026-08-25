## Usage: gdb -batch -x filename executable
catch throw
#r < run.txt
#r < trappe000_fstin.txt
#r < ljr001_fstin.txt
#r < patch001_fstin.txt
r --gtest_filter=BuildRecursiveTable*
bt
