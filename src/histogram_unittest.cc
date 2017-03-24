#include <gtest/gtest.h>
#include <limits.h>
#include "histogram.h"

// in order to test writing/reading restarts, and cloning,
//  this function was constructed to call 4 times
void testhist(Histogram *h) {
  EXPECT_NEAR(-2, h->min(), doubleTolerance);
  EXPECT_NEAR(6, h->max(), doubleTolerance);
  EXPECT_EQ(4, int(h->hist().size()));
  EXPECT_EQ(1, h->hist()[0]);
  EXPECT_EQ(1, h->hist()[1]);
  EXPECT_EQ(1, h->hist()[2]);
  EXPECT_EQ(1, h->hist()[3]);
  EXPECT_NEAR(-1, h->bin2m(0), doubleTolerance);
  EXPECT_NEAR(1, h->bin2m(1), doubleTolerance);
  EXPECT_NEAR(3, h->bin2m(2), doubleTolerance);
  EXPECT_NEAR(5, h->bin2m(3), doubleTolerance);
  EXPECT_NEAR(0, h->bin(-2), doubleTolerance);
  EXPECT_NEAR(0, h->bin(-doubleTolerance), doubleTolerance);
  EXPECT_NEAR(1, h->bin(0), doubleTolerance);
  EXPECT_NEAR(3, h->bin(5.5), doubleTolerance);
}

void testhistCenterZero(Histogram *h) {
  EXPECT_NEAR(-3, h->min(), doubleTolerance);
  EXPECT_NEAR(5, h->max(), doubleTolerance);
  EXPECT_EQ(4, int(h->hist().size()));
  EXPECT_EQ(1, h->hist()[0]);
  EXPECT_EQ(1, h->hist()[1]);
  EXPECT_EQ(1, h->hist()[2]);
  EXPECT_EQ(1, h->hist()[3]);
  EXPECT_NEAR(-2, h->bin2m(0), doubleTolerance);
  EXPECT_NEAR(0, h->bin2m(1), doubleTolerance);
  EXPECT_NEAR(2, h->bin2m(2), doubleTolerance);
  EXPECT_NEAR(4, h->bin2m(3), doubleTolerance);
  EXPECT_NEAR(0, h->bin(-2), doubleTolerance);
  EXPECT_NEAR(0, h->bin(-1-doubleTolerance), doubleTolerance);
  EXPECT_NEAR(1, h->bin(0), doubleTolerance);
  EXPECT_NEAR(3, h->bin(4), doubleTolerance);
}

TEST(Histogram, constructor) {
  { Histogram h(2);
    h.accumulate(3);
    EXPECT_NEAR(4, h.max(), doubleTolerance);
    EXPECT_NEAR(2, h.min(), doubleTolerance);
  }

  { Histogram h(2);
    h.accumulate(-3);
    EXPECT_NEAR(-2, h.max(), doubleTolerance);
    EXPECT_NEAR(-4, h.min(), doubleTolerance);
  }

  { Histogram h(2);
    h.accumulate(0);
    EXPECT_NEAR(2, h.max(), doubleTolerance);
    EXPECT_NEAR(0, h.min(), doubleTolerance);
    EXPECT_NEAR(0, h.bin(0), doubleTolerance);
    EXPECT_NEAR(1, h.bin2m(0), doubleTolerance);
    h.accumulate(5);
    EXPECT_NEAR(6, h.max(), doubleTolerance);
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
    EXPECT_NEAR(1, h.max(), doubleTolerance);
    EXPECT_NEAR(-1, h.min(), doubleTolerance);
    EXPECT_NEAR(0, h.bin(0), doubleTolerance);
    EXPECT_NEAR(0, h.bin2m(0), doubleTolerance);
    h.accumulate(4);
    EXPECT_NEAR(5, h.max(), doubleTolerance);
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

