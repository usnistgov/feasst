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
#include "pair_ideal.h"

TEST(PairIdeal, ideal) {
  feasst::Space s(3,0);
  const double boxl = 24.8586887;
  for (int dim=0; dim < s.dimen(); ++dim) s.lset(boxl,dim);
  const double rCut = 12.42934435;
  s.readXYZBulk(3, "water", "../unittest/spce/test52.xyz");
  feasst::PairIdeal p(&s, rCut);
  p.initEnergy();
  EXPECT_NEAR(0, p.peTot(), 1e-14);
}

