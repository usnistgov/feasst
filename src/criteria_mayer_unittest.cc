#include <gtest/gtest.h>
#include <limits.h>
#include "criteria_mayer.h"
#include "pair_squarewell.h"

TEST(CriteriaMayer, ref) {
  feasst::Space space(3, 0);
  space.lset(8);
  space.addMolInit("../forcefield/data.lj");
  vector<double> xAdd(3, 0.);
  space.xAdd = xAdd;
  space.addMol("../forcefield/data.lj");
  xAdd[0] = 1.1;
  space.xAdd = xAdd;
  space.addMol("../forcefield/data.lj");
  feasst::PairSquareWell pair(&space, 2.);
  pair.initData("../forcefield/data.lj");
  pair.rCutijset(0, 0, 2.);
  pair.initEnergy();

  const double beta = 1.;
  feasst::CriteriaMayer crit(beta);
  crit.initPairRef(&pair);
  crit.store(&space, &pair);
  EXPECT_EQ(0, crit.accept(0, 0., "move", 0));

}
