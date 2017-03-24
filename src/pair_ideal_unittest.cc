#include <gtest/gtest.h>
#include "pair_ideal.h"

TEST(PairIdeal, ideal) {
  Space s(3,0);
  const double boxl = 24.8586887;
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  const double rCut = 12.42934435;
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  PairIdeal p(&s, rCut);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 1e-14);
}

