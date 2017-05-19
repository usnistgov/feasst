#include <gtest/gtest.h>
#include <limits.h>
#include "criteria_mayer.h"
#include "pair_hs.h"

TEST(CriteriaMayer, ref) {
  feasst::Space space(3, 0);
  space.lset(8);
  space.addMolInit("../forcefield/data.lj");
  feasst::PairHS pair(&space, 1.);
  pair.initData("../forcefield/data.lj");

  const double beta = 1.;
  feasst::CriteriaMayer crit(beta);
  crit.initPairRef(&pair);
  EXPECT_EQ(0, crit.accept(0, 0., "move", 0));

}
