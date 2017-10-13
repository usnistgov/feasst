/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include <gtest/gtest.h>
#include <limits.h>
#include "random_nr3.h"

using namespace feasst;

TEST(Random, NR3ran) {
  RandomNR3 ran(17);
  EXPECT_EQ(269952321389814056u, ran.int64());
  EXPECT_EQ(7477734313819993120u, ran.int64());
  EXPECT_EQ(16294976781531816119u, ran.int64());
  EXPECT_EQ(17039904789424739738u, ran.int64());
  EXPECT_EQ(4945048831639962635u, ran.int64());
  EXPECT_EQ(1565409385732501729u, ran.int64());
  EXPECT_EQ(7095006703038622919u, ran.int64());
  EXPECT_EQ(13927236388846696772u, ran.int64());
  EXPECT_EQ(150171266583137103u, ran.int64());
  ran.writeRestart("tmp/rstnr3");
  RandomNR3 ran2("tmp/rstnr3");
  EXPECT_EQ(8092874815854888167u, ran.int64());
  EXPECT_EQ(8092874815854888167u, ran2.int64());
  for (int i = 0; i < 10; ++i) {
    EXPECT_EQ(ran.int64(), ran2.int64());
  }
}
