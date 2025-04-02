## Usage: gdb -batch -x filename executable
catch throw
#r < ljgibbs0_run.txt
r --gtest_filter=Configuration.type_to_file_name
bt
