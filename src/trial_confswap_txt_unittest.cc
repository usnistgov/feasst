#include <gtest/gtest.h>
#include "space.h"
#include "pair.h"
#include "pair_ideal.h"
#include "pair_lj.h"
#include "pair_patch_kf.h"
#include "pair_lj_coul_ewald.h"
#include "functions.h"
#include "criteria.h"
#include "criteria_metropolis.h"
#include "criteria_wltmmc.h"
#include "trial_transform.h"
#include "trial_confswap_txt.h"

using namespace feasst;

TEST(TrialConfSwapTXT, confswap) {
  Space s(3, 0);
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(12,dim);
  s.addMolInit("../forcefield/data.lj");
  for (int i = 0; i < 12; ++i) s.addMol("../forcefield/data.lj");
  PairLJ p(&s, 3);
  CriteriaMetropolis c(0.5, 0.01);
  TrialTransform tt(&s, &p, &c, "translate");
  tt.maxMoveParam = 5;
  TrialTransform tr(&s, &p, &c, "rotate");
  TrialConfSwapTXT tcs(&s, &p, &c);
  p.initEnergy();

  // clone
  shared_ptr<Space> s2 = s.cloneShrPtr();
  Pair* p2 = p.clone(s2.get());
  CriteriaMetropolis* c2 = c.clone();
  shared_ptr<Trial> tt2 = tt.cloneShrPtr(s2.get(), p2, c2);
  shared_ptr<Trial> tr2 = tr.cloneShrPtr(s2.get(), p2, c2);
  shared_ptr<TrialConfSwapTXT> tcs2 = tcs.cloneShrPtr(s2.get(), p2, c2);

  // initiate overlap
  tcs.initProc(0);
  tcs.initMType("nmol");
  tcs.addProcOverlap(s.nMol(), 1);
  tcs2->initProc(1);
  tcs2->initMType("nmol");
  tcs2->addProcOverlap(s.nMol(), 0);

  ranInitByDate();
  const int nAttempts = 300;
  for (int i = 0; i < nAttempts; ++i) {
    tt.attempt();
    tr.attempt();
    tcs.attempt();
    tt2->attempt();
    tr2->attempt();
    tcs2->attempt();
  }

  EXPECT_NE(tcs.attempted(), 0);
  EXPECT_NE(tcs2->attempted(), 0);
  EXPECT_NE(tcs.accepted(), 0);
  EXPECT_NE(tcs2->accepted(), 0);

  // free memory
  delete p2;
  delete c2;
}

