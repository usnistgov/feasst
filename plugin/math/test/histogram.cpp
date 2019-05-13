#include "utils/test/utils.h"
#include <limits.h>
#include "math/include/histogram.h"
#include "math/include/constants.h"
#include "utils/include/debug.h"

namespace feasst {

/// HWH fix this testing
/// HWH test setting bins by formula.

// in order to test writing/reading restarts, and cloning,
//  this function was constructed to call 4 times
void testhist(Histogram *h) {
  EXPECT_NEAR(-2, h->min(), NEAR_ZERO);
  EXPECT_NEAR(6, h->max(), NEAR_ZERO);
  EXPECT_EQ(4, int(h->histogram().size()));
  EXPECT_EQ(1, h->histogram()[0]);
  EXPECT_EQ(1, h->histogram()[1]);
  EXPECT_EQ(1, h->histogram()[2]);
  EXPECT_EQ(1, h->histogram()[3]);
  EXPECT_NEAR(-1, h->center_of_bin(0), NEAR_ZERO);
  EXPECT_NEAR(1, h->center_of_bin(1), NEAR_ZERO);
  EXPECT_NEAR(3, h->center_of_bin(2), NEAR_ZERO);
  EXPECT_NEAR(5, h->center_of_bin(3), NEAR_ZERO);
  EXPECT_NEAR(0, h->bin(-2), NEAR_ZERO);
  EXPECT_NEAR(0, h->bin(-1e-15), NEAR_ZERO);
  EXPECT_NEAR(1, h->bin(0), NEAR_ZERO);
  EXPECT_NEAR(3, h->bin(5.5), NEAR_ZERO);
}

void testhistCenterZero(Histogram *h) {
  EXPECT_NEAR(-3, h->min(), NEAR_ZERO);
  EXPECT_NEAR(5, h->max(), NEAR_ZERO);
  EXPECT_EQ(4, int(h->histogram().size()));
  EXPECT_EQ(1, h->histogram()[0]);
  EXPECT_EQ(1, h->histogram()[1]);
  EXPECT_EQ(1, h->histogram()[2]);
  EXPECT_EQ(1, h->histogram()[3]);
  EXPECT_NEAR(-2, h->center_of_bin(0), NEAR_ZERO);
  EXPECT_NEAR(0, h->center_of_bin(1), NEAR_ZERO);
  EXPECT_NEAR(2, h->center_of_bin(2), NEAR_ZERO);
  EXPECT_NEAR(4, h->center_of_bin(3), NEAR_ZERO);
  EXPECT_NEAR(0, h->bin(-2), NEAR_ZERO);
  EXPECT_NEAR(0, h->bin(-1-1e-15), NEAR_ZERO);
  EXPECT_NEAR(1, h->bin(0), NEAR_ZERO);
  EXPECT_NEAR(3, h->bin(4), NEAR_ZERO);
}

TEST(Histogram, constructor) {
  { Histogram h;
    h.set_width_center(2, 3);
    h.add(3);
    EXPECT_NEAR(4, h.max(), NEAR_ZERO);
    EXPECT_NEAR(2, h.min(), NEAR_ZERO);
  }

  { Histogram h;
    h.set_width_center(2, -3);
    h.add(-3);
    EXPECT_NEAR(-2, h.max(), NEAR_ZERO);
    EXPECT_NEAR(-4, h.min(), NEAR_ZERO);
  }

  { Histogram h;
    h.set_width_center(2, 1);
    h.add(0);
    EXPECT_NEAR(2, h.max(), NEAR_ZERO);
    EXPECT_NEAR(0, h.min(), NEAR_ZERO);
    EXPECT_NEAR(0, h.bin(0), NEAR_ZERO);
    EXPECT_NEAR(1, h.center_of_bin(0), NEAR_ZERO);
    h.add(5);
    EXPECT_NEAR(6, h.max(), NEAR_ZERO);
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

  { Histogram h;
    h.set_width_center(2, 0);
    h.add(0);
    EXPECT_NEAR(1, h.max(), NEAR_ZERO);
    EXPECT_NEAR(-1, h.min(), NEAR_ZERO);
    EXPECT_NEAR(0, h.bin(0), NEAR_ZERO);
    EXPECT_NEAR(0, h.center_of_bin(0), NEAR_ZERO);
    h.add(4);
    EXPECT_NEAR(5, h.max(), NEAR_ZERO);
    EXPECT_EQ(3, int(h.size()));
    EXPECT_EQ(1, h.histogram()[0]);
    EXPECT_EQ(0, h.histogram()[1]);
    EXPECT_EQ(1, h.histogram()[2]);

    h.add(2);
    EXPECT_EQ(1, h.histogram()[1]);

    h.add(-2);

    testhistCenterZero(&h);

    // serialize
    Histogram h2 = test_serialize(h);
    testhistCenterZero(&h2);
    EXPECT_EQ(1, h2.histogram()[1]);
  }
}

}  // namespace feasst
