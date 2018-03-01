/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include "arguments.h"
#include "functions.h"

TEST(Arguments, Arguments) {
  feasst::Arguments args;
  args.initArgs("test", {{"key1", "val1"}});
  EXPECT_EQ(args.key("yo").dflt("hi").str(), "hi");

  // test if a key is empty before it is set
  try {
    args.empty();
    CATCH_PHRASE("key must be set before");
  }

  // test if a key is empty after set
  EXPECT_TRUE(args.key("yo").empty());

  // test if all keywords have been used
  try {
    EXPECT_FALSE(args.checkAllArgsUsed());
    CATCH_PHRASE("All keywords provided in args must be used");
  }

  // test if provided key is empty
  EXPECT_FALSE(args.key("key1").empty());

  // check second key pair without providing key, because already set
  EXPECT_EQ(args.dflt("hi").str(), "val1");

  // alternatively, set key and check simultaneously
  EXPECT_EQ(args.key("key1").dflt("hi").str(), "val1");

  // test if key1 is removed using rm()
  const int nargs = args.size();
  EXPECT_EQ(args.key("key1").dflt("hi").rm().str(), "val1");
  EXPECT_EQ(args.size(), nargs - 1);

  // check if all args were used
  EXPECT_TRUE(args.checkAllArgsUsed());
}

TEST(Arguments, integer) {
  try {
    feasst::Arguments args;
    args.initArgs("test", {{"key1", "val1"}});
    args.key("key1").integer();
    CATCH_PHRASE("was expected to be an integer");
  }
  try {
    feasst::Arguments args;
    args.initArgs("test", {{"key1", "1.1"}});
    args.key("key1").integer();
    CATCH_PHRASE("was expected to be an integer");
  }

  feasst::Arguments args;
  args.initArgs("test", {{"key1", "1."}});
  args.key("key1").integer();
}

TEST(Arguments, dble) {
  try {
    feasst::Arguments args;
    args.initArgs("test", {{"key1", "mymypie"}});
    args.key("key1").dble();
    CATCH_PHRASE("was expected to be a double precision floating point number");
  }
  feasst::Arguments args;
  args.initArgs("test", {{"key1", "3.1415"}});
  EXPECT_NEAR(3.1415, args.key("key1").dble(), feasst::DTOL);
}

TEST(Arguments, arglist) {
  feasst::arglist argls = {{ {"set1", {{"key1", "val1"}} },
                             {"set2", {{"key2", "val2"}} } }};
  auto set = argls.find("set2");
  auto pair = set->second.find("key2");
  EXPECT_NE(set->second.end(), pair);
  EXPECT_TRUE(pair->second == "val2");
  cout << pair->second << endl;
}
