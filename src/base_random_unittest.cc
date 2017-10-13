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
#include "base_random.h"
#include "accumulator.h"
#include "histogram.h"

using namespace feasst;

TEST(BaseRandom, uniformRanNum) {
  ranInitByDate();
  BaseRandom ran;
  const int min = 5, max = 8, n = 10000;
  vector<int> x(max);
  for (int i = 0; i < n; ++i) {
    ++x.at(ran.uniformRanNum(min, max) - 1);
  }
  for (int i = min; i <= max; ++i) {
    EXPECT_NEAR(1./(max - min + 1), double(x.at(i-1))/n, 5e-2);
  }
  EXPECT_EQ(0, x.at(min-2));
  ran.writeRngRestart("tmp/rstnr3");
  BaseRandom ran2("tmp/rstnr3");
  EXPECT_NEAR(ran.uniformRanNum(), ran2.uniformRanNum(), 1e-19);
}

TEST(BaseRandom, stdNormRanNum) {
  ranInitByDate();
  BaseRandom ran;
  const int n = 10000;
  Accumulator a;
  for (int i = 0; i < n; ++i) {
    const double random = ran.stdNormRanNum();
    a.accumulate(random);
  }
  EXPECT_NEAR(0, a.average(), 5e-2);
  EXPECT_NEAR(1, a.stdev(), 5e-2);

}

TEST(BaseRandom, gaussRanNum) {
  ranInitByDate();
  BaseRandom ran;
  const int n = 10000;
  const double sig = 5, av = 10;
  Accumulator a;
  Histogram hist(0.5);
  for (int i = 0; i < n; ++i) {
    const double random = ran.gaussRanNum(sig, av);
    a.accumulate(random);
    hist.accumulate(random);
  }
  EXPECT_NEAR(av, a.average(), sig*5e-2);
  EXPECT_NEAR(sig, a.stdev(), sig*5e-2);
  hist.write("tmp/histhist");
}

TEST(BaseRandom, hash) {
  BaseRandom ran;
  std::string hash = ran.randomHash();
  EXPECT_NE(hash, ran.randomHash());
  EXPECT_NE(hash, ran.randomHash());
}
