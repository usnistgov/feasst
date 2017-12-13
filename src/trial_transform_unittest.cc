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

TEST(TrialTransform, args) {
  feasst::Space space;
  feasst::PairHardSphere pair(&space);
  feasst::CriteriaMetropolis crit(1., 1.);

  {
    feasst::TrialTransform trial(&pair, &crit,
      {{"type", "translate"},
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
       {"type", "meh"}});
    CATCH_PHRASE("(meh) not recognized");
  }

  try {
    feasst::TrialTransform trial(&pair, &crit,
      {{"/not/a/param", "error"},
       {"type", "translate"}});
    CATCH_PHRASE("is not recognized");
  }
}


