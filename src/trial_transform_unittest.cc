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
#include "pair_hard_sphere.h"
#include "criteria_metropolis.h"
#include "trial_transform.h"
#include "./group.h"

TEST(TrialTransform, args) {
  feasst::Space space;
  feasst::PairHardSphere pair(&space);
  feasst::CriteriaMetropolis crit(1., 1.);

  {
    feasst::TrialTransform trial(&pair, &crit,
      {{"transType", "translate"},
       {"maxMoveParam", "0.1"}});
    EXPECT_TRUE(trial.transType() == "translate");
    EXPECT_TRUE(trial.maxMoveParam == 0.1);
  }

  try {
    feasst::TrialTransform trial(&pair, &crit, {{"/not/a/param", "error"}});
    CATCH_PHRASE("is required for args");
    //CATCH_PHRASE("is not recognized");
  }

  try {
    feasst::TrialTransform trial(&pair, &crit,
      {{"/not/a/param", "error"},
       {"transType", "meh"}});
    CATCH_PHRASE("(meh) not recognized");
  }

  try {
    feasst::TrialTransform trial(&pair, &crit,
      {{"/not/a/param", "error"},
       {"transType", "translate"}});
    CATCH_PHRASE("is not recognized");
  }
}

TEST(TrialTransform, group) {
  feasst::Space space;
  auto g2 = make_shared<feasst::Group>();
  g2->molid(1);
  g2->initName("lj");
  space.initGroup(g2);
  feasst::PairHardSphere pair(&space);
  pair.initData("../forcefield/data.atom");
  pair.initData("../forcefield/data.lj");
  pair.addMol("../forcefield/data.atom");
  pair.addMol("../forcefield/data.lj");
  feasst::CriteriaMetropolis crit(1., 1.);
  feasst::TrialTransform trial(&pair, &crit,
    {{"transType", "translate"},
     {"maxMoveParam", "0.1"}});
  trial.initGroup("lj");
  pair.printXYZ("tmp/hitt", 1);
  for (int i = 0; i < 12; ++i) {
    trial.attempt();
    pair.printXYZ("tmp/hitt", 0);
  }
}

