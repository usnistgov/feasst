# This script performs a build and full gambit of tests
# 1. compile
# 2. valgrind
# 3. C++ unittests
# 4. python unittests
# 5. coverage
# 6. documentation
# 7. test cases
# compilation with cmake and make
../buildexample/compile.sh
../buildexample/unittest.sh   # this can't run at the same time as valgrind
../buildexample/makecov.sh    # this can't run at the same time as valgrind

# run more tests in parallel
pidArr=()
for tests in valgrind swigtest cpplint pylint makehtml; do
  ../buildexample/$tests.sh &
  pidArr+=($!)
done
wait ${pidArr[@]}

echo "Nightly Build completed "`date` >> summary.log

