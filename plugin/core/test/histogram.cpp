#include <gtest/gtest.h>
#include <limits.h>
#include "core/include/histogram.h"

/// HWH test setting bins by formula.

// in order to test writing/reading restarts, and cloning,
//  this function was constructed to call 4 times
void testhist(feasst::Histogram *h) {
  EXPECT_NEAR(-2, h->min(), 1e-15);
  EXPECT_NEAR(6, h->max(), 1e-15);
  EXPECT_EQ(4, int(h->histogram().size()));
  EXPECT_EQ(1, h->histogram()[0]);
  EXPECT_EQ(1, h->histogram()[1]);
  EXPECT_EQ(1, h->histogram()[2]);
  EXPECT_EQ(1, h->histogram()[3]);
  EXPECT_NEAR(-1, h->center_of_bin(0), 1e-15);
  EXPECT_NEAR(1, h->center_of_bin(1), 1e-15);
  EXPECT_NEAR(3, h->center_of_bin(2), 1e-15);
  EXPECT_NEAR(5, h->center_of_bin(3), 1e-15);
  EXPECT_NEAR(0, h->bin(-2), 1e-15);
  EXPECT_NEAR(0, h->bin(-1e-15), 1e-15);
  EXPECT_NEAR(1, h->bin(0), 1e-15);
  EXPECT_NEAR(3, h->bin(5.5), 1e-15);
}

void testhistCenterZero(feasst::Histogram *h) {
  EXPECT_NEAR(-3, h->min(), 1e-15);
  EXPECT_NEAR(5, h->max(), 1e-15);
  EXPECT_EQ(4, int(h->histogram().size()));
  EXPECT_EQ(1, h->histogram()[0]);
  EXPECT_EQ(1, h->histogram()[1]);
  EXPECT_EQ(1, h->histogram()[2]);
  EXPECT_EQ(1, h->histogram()[3]);
  EXPECT_NEAR(-2, h->center_of_bin(0), 1e-15);
  EXPECT_NEAR(0, h->center_of_bin(1), 1e-15);
  EXPECT_NEAR(2, h->center_of_bin(2), 1e-15);
  EXPECT_NEAR(4, h->center_of_bin(3), 1e-15);
  EXPECT_NEAR(0, h->bin(-2), 1e-15);
  EXPECT_NEAR(0, h->bin(-1-1e-15), 1e-15);
  EXPECT_NEAR(1, h->bin(0), 1e-15);
  EXPECT_NEAR(3, h->bin(4), 1e-15);
}

TEST(Histogram, constructor) {
  { feasst::Histogram h;
    h.set_width_center(2, 3);
    h.add(3);
    EXPECT_NEAR(4, h.max(), 1e-15);
    EXPECT_NEAR(2, h.min(), 1e-15);
  }

  { feasst::Histogram h;
    h.set_width_center(2, -3);
    h.add(-3);
    EXPECT_NEAR(-2, h.max(), 1e-15);
    EXPECT_NEAR(-4, h.min(), 1e-15);
  }

  { feasst::Histogram h;
    h.set_width_center(2, 1);
    h.add(0);
    EXPECT_NEAR(2, h.max(), 1e-15);
    EXPECT_NEAR(0, h.min(), 1e-15);
    EXPECT_NEAR(0, h.bin(0), 1e-15);
    EXPECT_NEAR(1, h.center_of_bin(0), 1e-15);
    h.add(5);
    EXPECT_NEAR(6, h.max(), 1e-15);
    EXPECT_EQ(3, int(h.size()));
    EXPECT_EQ(1, h.histogram()[0]);
    EXPECT_EQ(0, h.histogram()[1]);
    EXPECT_EQ(1, h.histogram()[2]);

    h.add(3);
    EXPECT_EQ(1, h.histogram()[1]);

    h.add(-1);

    testhist(&h);

//    h.writeRestart("tmp/tmph");
//    Histogram h2("tmp/tmph");
//    Histogram* h3 = h.clone();
//    shared_ptr<Histogram> h4 = h.cloneShrPtr();
//
//    testhist(&h2);
//    testhist(h3);
//    testhist(h4.get());
//
//    delete h3;
  }

  { feasst::Histogram h;
    h.set_width_center(2, 0);
    h.add(0);
    EXPECT_NEAR(1, h.max(), 1e-15);
    EXPECT_NEAR(-1, h.min(), 1e-15);
    EXPECT_NEAR(0, h.bin(0), 1e-15);
    EXPECT_NEAR(0, h.center_of_bin(0), 1e-15);
    h.add(4);
    EXPECT_NEAR(5, h.max(), 1e-15);
    EXPECT_EQ(3, int(h.size()));
    EXPECT_EQ(1, h.histogram()[0]);
    EXPECT_EQ(0, h.histogram()[1]);
    EXPECT_EQ(1, h.histogram()[2]);

    h.add(2);
    EXPECT_EQ(1, h.histogram()[1]);

    h.add(-2);

    testhistCenterZero(&h);

//    h.writeRestart("tmp/tmph");
//    Histogram h2("tmp/tmph");
//    Histogram* h3 = h.clone();
//    shared_ptr<Histogram> h4 = h.cloneShrPtr();
//
//    testhistCenterZero(&h2);
//    testhistCenterZero(h3);
//    testhistCenterZero(h4.get());
//
//    delete h3;
  }
}

