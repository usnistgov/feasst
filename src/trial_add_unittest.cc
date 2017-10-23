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
#include "pair_lj_multi.h"
#include "trial_transform.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "mc.h"

using namespace feasst;

TEST(TrialAdd, confine) {
  Space space(3, 0);  // create 3D space
  space.lset(10);     // create cubic PBC box length of 10, origin at center
  stringstream addMol;
  addMol << space.install_dir() << "/forcefield/data.lj";
  stringstream addWall;
  addWall << space.install_dir() << "/forcefield/data.ljb";

  // Initialize LJ interactions
  const double rCut = 3.;
  PairLJMulti pair(&space, rCut);
  pair.initData(addMol.str());
  pair.initData(addWall.str());
  pair.epsijset(1, 1, 0.);  //!< turn off Wall-Wall interactions
  pair.linearShift();

  // create wall of frozen LJ particles at z=0 plane
  vector<double> xAdd(space.dimen());
  const double latticeSpacing = 0.5;
  for (double x = -0.5*space.l(0); x < 0.5*space.l(0) - 0.0001;
       x += latticeSpacing) {
    xAdd[0] = x;
    for (double y = -0.5*space.l(1); y < 0.5*space.l(1) - 0.0001;
         y += latticeSpacing) {
      xAdd[1] = y;
      space.xAdd = xAdd;
      space.addMol(addWall.str());
    }
  }
  pair.addPart();  //!< tell pair class that you updated particles

  // Initialize GCMC simulation
  pair.initEnergy();
  const double beta = 1., activ = exp(-7.);
  CriteriaMetropolis crit(beta, activ);
  MC mc(&space, &pair, &crit);

  // Initialize trial moves
  shared_ptr<TrialTransform> ttrans = make_shared<TrialTransform>("translate");
  ttrans->selectType(addMol.str().c_str());
  mc.weight = 1.;
  mc.initTrial(ttrans);

  shared_ptr<TrialDelete> tdel = make_shared<TrialDelete>(addMol.str().c_str());
  mc.weight = 0.5;
  mc.initTrial(tdel);

  shared_ptr<TrialAdd> tadd = make_shared<TrialAdd>(addMol.str().c_str());
  mc.weight = 0.5;
  mc.initTrial(tadd);

  // Note that confine can't be called until the space object pointer is
  // initialized. This means it must follow "initTrial" in this example.
  const double upper = 2., lower = 0.;
  const int confineDim = 2;  // z-dimension
  tdel->confine(upper, lower, confineDim);
  tadd->confine(upper, lower, confineDim);

  // Run MC simulation
  const int nPrint = 1e1;
  mc.initLog("tmp/ljaddconfine", nPrint);
  mc.initMovie("tmp/ljaddconfine", nPrint);
  mc.initRestart("tmp/ljaddconfinerst", nPrint);
  mc.setNFreqCheckE(nPrint, 1e-4);
  mc.setNFreqTune(nPrint);
  mc.runNumTrials(1e2);
}

