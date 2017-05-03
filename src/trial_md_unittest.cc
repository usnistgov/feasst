#include <gtest/gtest.h>
#include "trial_md.h"
#include "pair_lj_multi.h"
#include "criteria_metropolis.h"

using namespace feasst;

TEST(TrialMD, velocity) {
  Space space(3, 0);
  space.addMolInit("../forcefield/data.atom");
  for (int i = 0; i < 10; ++i) {
    space.addMol("../forcefield/data.atom");
  }
  PairLJMulti pair(&space, 3.);
  CriteriaMetropolis criteria(0.5, 0.);
  TrialMD tmd(&space, &pair, &criteria);

  tmd.initMomentum();
  EXPECT_NEAR(tmd.temperature(), 1./criteria.beta(), 100*doubleTolerance);
  //tmd.attempt();
}


