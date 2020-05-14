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
#include "pair_ideal.h"
#include "pair_hard_sphere.h"
#include "pair_lj.h"
#include "pair_lj_coul_ewald.h"
#include "mc_wltmmc.h"
#include "ui_abbreviated.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_transform.h"

TEST(WLTMMC, windowing) {
  feasst::ranInitByDate();
  auto space = feasst::makeSpace(
    {{"dimen", "3"},
     {"boxLength", "8"}});
  auto pair = feasst::makePairLJ(space,
    {{"rCut", "3"},
     {"cutType", "lrc"},
     {"molTypeInForcefield", "data.lj"}});
  auto criteria = feasst::makeCriteriaWLTMMC(
    {{"beta", "1"},
     {"activ", feasst::str(exp(-1))},
     {"mType", "nmol"},
     {"nMax", "60"},
     {"nMin", "10"}});
  feasst::WLTMMC mc(pair, criteria);
  try {
    mc.initWindows({-10, 70});
    CATCH_PHRASE("midpoint -10 is less than minimum");
  }
  try {
    mc.initWindows({20, 70});
    CATCH_PHRASE("midpoint 70 is greater than maximum");
  }
  try {
    mc.initWindows({20, 50.5});
    CATCH_PHRASE("midPoint 50.5 does not lie in the middle");
  }
}
