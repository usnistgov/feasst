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
  TrialMD tmd(&pair, &criteria);

  tmd.initMomentum();
  EXPECT_NEAR(tmd.temperature(), 1./criteria.beta(), 100*DTOL);
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
  EXPECT_NEAR(utot, p.peTot() + tmd->kineticEnergy(), 0.75);
}


