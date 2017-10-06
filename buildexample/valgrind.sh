file=$(basename $0 .sh)
valgrind --suppressions=../tools/valgrind/gomp.supp --suppressions=../tools/valgrind/gomp2.supp --suppressions=../tools/valgrind/gomp3.supp --suppressions=../tools/valgrind/gomp4.supp --suppressions=../tools/valgrind/gomp5.supp --track-origins=yes --leak-check=full --show-leak-kinds=all ./bin/unittest --gtest_shuffle > $file.log 2>&1
echo "**** $file ****" >> summary.log
tail -4 $file.log >> summary.log
echo "" >> summary.log
