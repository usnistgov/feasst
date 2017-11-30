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

  // test if all keywords have been used
  try {
    EXPECT_FALSE(args.checkAllArgsUsed());
    CATCH_PHRASE("All keywords provided in args must be used");
  }

  // use keyword before args destructor to prevent exception
  EXPECT_EQ(args.key("key1").dflt("hi").str(), "val1");
  EXPECT_TRUE(args.checkAllArgsUsed());
}

