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
  ../tools/build/${task}.sh
done

# run more tests in parallel
pidArr=()
for tests in cpplint pylint makehtml; do
  ../tools/build/$tests.sh &
  pidArr+=($!)
done
wait ${pidArr[@]}

echo "Nightly Build completed "`date` >> summary.log

