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
#include "histogram.h"

using namespace feasst;

// in order to test writing/reading restarts, and cloning,
//  this function was constructed to call 4 times
void testhist(Histogram *h) {
  EXPECT_NEAR(-2, h->min(), DTOL);
  EXPECT_NEAR(6, h->max(), DTOL);
  EXPECT_EQ(4, int(h->hist().size()));
  EXPECT_EQ(1, h->hist()[0]);
  EXPECT_EQ(1, h->hist()[1]);
  EXPECT_EQ(1, h->hist()[2]);
  EXPECT_EQ(1, h->hist()[3]);
  EXPECT_NEAR(-1, h->bin2m(0), DTOL);
  EXPECT_NEAR(1, h->bin2m(1), DTOL);
  EXPECT_NEAR(3, h->bin2m(2), DTOL);
  EXPECT_NEAR(5, h->bin2m(3), DTOL);
  EXPECT_NEAR(0, h->bin(-2), DTOL);
  EXPECT_NEAR(0, h->bin(-DTOL), DTOL);
  EXPECT_NEAR(1, h->bin(0), DTOL);
  EXPECT_NEAR(3, h->bin(5.5), DTOL);
}

void testhistCenterZero(Histogram *h) {
  EXPECT_NEAR(-3, h->min(), DTOL);
  EXPECT_NEAR(5, h->max(), DTOL);
  EXPECT_EQ(4, int(h->hist().size()));
  EXPECT_EQ(1, h->hist()[0]);
  EXPECT_EQ(1, h->hist()[1]);
  EXPECT_EQ(1, h->hist()[2]);
  EXPECT_EQ(1, h->hist()[3]);
  EXPECT_NEAR(-2, h->bin2m(0), DTOL);
  EXPECT_NEAR(0, h->bin2m(1), DTOL);
  EXPECT_NEAR(2, h->bin2m(2), DTOL);
  EXPECT_NEAR(4, h->bin2m(3), DTOL);
  EXPECT_NEAR(0, h->bin(-2), DTOL);
  EXPECT_NEAR(0, h->bin(-1-DTOL), DTOL);
  EXPECT_NEAR(1, h->bin(0), DTOL);
  EXPECT_NEAR(3, h->bin(4), DTOL);
}

TEST(Histogram, constructor) {
  { Histogram h(2);
    h.accumulate(3);
    EXPECT_NEAR(4, h.max(), DTOL);
    EXPECT_NEAR(2, h.min(), DTOL);
  }

  { Histogram h(2);
    h.accumulate(-3);
    EXPECT_NEAR(-2, h.max(), DTOL);
    EXPECT_NEAR(-4, h.min(), DTOL);
  }

  { Histogram h(2);
    h.accumulate(0);
    EXPECT_NEAR(2, h.max(), DTOL);
    EXPECT_NEAR(0, h.min(), DTOL);
    EXPECT_NEAR(0, h.bin(0), DTOL);
    EXPECT_NEAR(1, h.bin2m(0), DTOL);
    h.accumulate(5);
    EXPECT_NEAR(6, h.max(), DTOL);
    EXPECT_EQ(3, int(h.hist().size()));
    EXPECT_EQ(1, h.hist()[0]);
    EXPECT_EQ(0, h.hist()[1]);
    EXPECT_EQ(1, h.hist()[2]);

    h.accumulate(3);
    EXPECT_EQ(1, h.hist()[1]);

    h.accumulate(-1);

    testhist(&h);

    h.writeRestart("tmp/tmph");
    Histogram h2("tmp/tmph");
    Histogram* h3 = h.clone();
    shared_ptr<Histogram> h4 = h.cloneShrPtr();

    testhist(&h2);
    testhist(h3);
    testhist(h4.get());

    delete h3;
  }

  { Histogram h(2);
    h.centerZero();
    h.accumulate(0);
    EXPECT_NEAR(1, h.max(), DTOL);
    EXPECT_NEAR(-1, h.min(), DTOL);
    EXPECT_NEAR(0, h.bin(0), DTOL);
    EXPECT_NEAR(0, h.bin2m(0), DTOL);
    h.accumulate(4);
    EXPECT_NEAR(5, h.max(), DTOL);
    EXPECT_EQ(3, int(h.hist().size()));
    EXPECT_EQ(1, h.hist()[0]);
    EXPECT_EQ(0, h.hist()[1]);
    EXPECT_EQ(1, h.hist()[2]);

    h.accumulate(2);
    EXPECT_EQ(1, h.hist()[1]);

    h.accumulate(-2);

    testhistCenterZero(&h);

    h.writeRestart("tmp/tmph");
    Histogram h2("tmp/tmph");
    Histogram* h3 = h.clone();
    shared_ptr<Histogram> h4 = h.cloneShrPtr();

    testhistCenterZero(&h2);
    testhistCenterZero(h3);
    testhistCenterZero(h4.get());

    delete h3;
  }
}

