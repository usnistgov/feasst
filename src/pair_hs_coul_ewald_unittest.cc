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
#include "pair_hs_coul_ewald.h"

TEST(PairHSCoulEwald, RPM) {
  auto space = feasst::makeSpace(
   {{"dimen", "3"},
    {"boxLength", "20"}});

  auto pair = feasst::makePairHSCoulEwald(space.get(), {{"rCut", "10"}});
  pair->initData("../forcefield/data.rpm_plus");
  pair->initData("../forcefield/data.rpm_minus");
  pair->initKSpace(5.6, 38);
  for (int iType = 0; iType < space->nParticleTypes(); ++iType) {
    for (int jType = 0; jType < space->nParticleTypes(); ++jType) {
      pair->rCutijset(iType, jType, pair->rCut());
    }
  }

  // test hard sphere
  vector<double> xAdd(space->dimen(), 0.);  // vector position at origin
  pair->addMol(xAdd, "../forcefield/data.rpm_plus");
  pair->addPart();  // update kspace
  xAdd[0] = 0.99;   // set x-coordinate of next particle
  pair->addMol(xAdd, "../forcefield/data.rpm_minus");
  pair->addPart();
  pair->initEnergy(); // calculate energy
  pair->printXYZ("hi", 1);  // remove this eventually
  EXPECT_GT(pair->peTot(), 1e200);
  EXPECT_NEAR(0, pair->peLRC(), feasst::DTOL);

  // test charge
  xAdd[0] = 1.01 - space->x(1, 0);  // preparing to move particle to x=1.01
  space->transMol(1,  // translate mol 1 (second mol)
    xAdd);  // by vector displacement in xAdd
  pair->initEnergy();
  pair->printXYZ("hi", 0);
  EXPECT_LT(pair->peTot(), 1e200);
  EXPECT_NEAR(0, pair->peLRC(), feasst::DTOL);

  // calculate the charge analytically by expanding the box, 5.5.2 of allen and tildesley
}

