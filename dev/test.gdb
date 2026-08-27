## Usage: export DEBUGINFOD_URLS=""; gdb -batch -x filename executable
catch throw
#r < chemi000_fstin.txt
r < runs/temp373/1w/mu_1000/acid_zeo000_fstin.txt
#--gtest_filter=BuildRecursiveTable*
bt
