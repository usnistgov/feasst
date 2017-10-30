# This script performs a build and full gambit of tests
# 1. compile
# 2. valgrind
# 3. C++ unittests
# 4. python unittests
# 5. coverage
# 6. documentation
# 7. test cases
# compilation with cmake and make

for task in compile unittest makecov valgrind testcase testcasepy; do
  ../buildexample/${task}.sh
done

# run more tests in parallel
pidArr=()
for tests in swigtest cpplint pylint makehtml; do
  ../buildexample/$tests.sh &
  pidArr+=($!)
done
wait ${pidArr[@]}

echo "Nightly Build completed "`date` >> summary.log

