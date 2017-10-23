#!/usr/bin/env bash
# converts C++ main file to gtest file
# usage ./cc2gtest [fileName]

SRC=$1
if [ -z $2 ]; then
  DES="$(basename $SRC .cc)"_unittest.cc
else
  DES=$2
fi

# write a macro which calls GTEST macro
cat << EOF > $DES
# define ASSERT(condition, message) \
if (! (condition)) { \
  std::stringstream err_msg; \
  std::cout << "# Assertion `" #condition "` failed in " << __FILE__ \
              << " line " << __LINE__ << ": " << message << std::endl; \
} \
EXPECT_TRUE(condition);
#include <gtest/gtest.h>
EOF

cat $SRC >> $DES

# remove getopt
awk 'f&&/GETOPT/{print "{";f=0} !f; /parse command-line arguments using/{f=1}' $DES > /tmp/tmp
mv /tmp/tmp $DES

# to begin, replace int main with GTEST macro
sed 's/^int main(\(.*\)) {  \/\/ \(.*\)/TEST(\2) {/' $DES > /tmp/tmp; mv /tmp/tmp $DES

# now, replace ASSERT with GTEST macro
sed 's/ASSERT(/CUSTOM_ASSERT(/' $DES > /tmp/tmp; mv /tmp/tmp $DES

