#include <gtest/gtest.h>
#include "trial_md.h"
#include "pair_lj_multi.h"
#include "criteria_metropolis.h"
#include "mc.h"

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

TEST(MCPrivate, ljMDconserve) {
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(6, dim);
  stringstream addMolType;
  addMolType << "../forcefield/data.lj";
  s.addMolInit(addMolType.str().c_str());
  PairLJ p(&s, 3.);
  p.initEnergy();
  CriteriaMetropolis c(1/0.85, 0.5);
  MC mc(&s, &p, &c);
  mc.nMolSeek(50, addMolType.str().c_str());
  shared_ptr<TrialMD> tmd = make_shared<TrialMD>();
  mc.initTrial(tmd);
  tmd->rescaleTemp = 0;
  mc.runNumTrials(1);
  const double utot = p.peTot() + tmd->kineticEnergy();
  mc.runNumTrials(200);
  EXPECT_NEAR(utot, p.peTot() + tmd->kineticEnergy(), 0.55);
}


