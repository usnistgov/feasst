## Usage: export DEBUGINFOD_URLS=""; gdb -batch -x filename executable
catch throw
#r < chemi000_fstin.txt
r --gtest_filter=MonteCarlo.rigid_body_connector
bt
