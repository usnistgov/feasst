/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include <gtest/gtest.h>
#include "space.h"
#include "pair.h"
#include "pair_ideal.h"
#include "pair_lj.h"
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
  s.initBoxLength(12);
  PairLJ p(&s, 3, {{"molType", "../forcefield/data.lj"}});
  for (int i = 0; i < 12; ++i) {
    p.addMol("../forcefield/data.lj");
  }
  CriteriaMetropolis c(0.5, 0.01);
  TrialTransform tt(&p, &c, "translate");
  tt.maxMoveParam = 5;
  TrialTransform tr(&p, &c, "rotate");
  TrialConfSwapTXT tcs(&p, &c);
  p.initEnergy();

  // clone
  shared_ptr<Space> s2 = s.cloneShrPtr();
  Pair* p2 = p.clone(s2.get());
  CriteriaMetropolis* c2 = c.clone();
  shared_ptr<Trial> tt2 = tt.cloneShrPtr(p2, c2);
  shared_ptr<Trial> tr2 = tr.cloneShrPtr(p2, c2);
  shared_ptr<TrialConfSwapTXT> tcs2 = tcs.cloneShrPtr(p2, c2);

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

