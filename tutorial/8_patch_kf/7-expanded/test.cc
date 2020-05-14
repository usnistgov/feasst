/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "feasst.h"

int main() {  // PatchKF, expandedNVT
  const int nMolMax = 48;
  const double betaMin = 10, betaMax = 100;
  const int betaBins = 11;
  feasst::ranInitByDate();
  auto space = feasst::makeSpace(
    {{"dimen", "2"},
     {"boxLength", "40."}});
  auto pair = feasst::makePairPatchKFMulti(space.get(),
    {{"rCut", "1.5"},  // Note: rCut is the patch well center-center distance
     {"patchAngle", "5"}});
  // Note: data file sets the bond length, epsilon (depth) and sigma
  // Pair Coeffs are in order [type] [eps] [sig]
  // Patch unit vectors have 0 sigma
  // Center of mass rotation centers of 0 eps (only in cg2_2patch)
  std::stringstream ssA;
  ssA << space->install_dir() << "/forcefield/data.cg4_3patch";
  pair->initData(ssA.str());
  std::stringstream ssB;
  //ssB << space->install_dir() << "/forcefield/data.cg3_2patch_linear";
  ssB << space->install_dir() << "/forcefield/data.cg2_2patch_linear";
  //ssB << space->install_dir() << "/forcefield/data.cg1_2patch_linear";
  pair->initData(ssB.str());

  // turn off the patchy interactions between like molecules
  pair->epsijset(2, 2, 0.);
  pair->epsijset(5, 5, 0.);

  // If you change the sigmas, rCut is not set dependent on sigma, so this
  // initIJ may need to be changed if you want rcut-sigma=constant
  pair->initIJ();

  space->updateCells(pair->rCut());

  auto criteria = feasst::makeCriteriaWLTMMC(
    {{"beta", feasst::str(betaMax)},   // beta=1/kT
     {"activ", feasst::str(exp(-7))},  // lnz = beta*mu - 3ln(debroglie)
     {"mType", "beta"},
     {"mMinCenter", feasst::str(betaMin)},
     {"mMaxCenter", feasst::str(betaMax)},
     {"nBin", feasst::str(betaBins)}
    });
  criteria->addActivity(exp(-7));
  criteria->collectInit();
  criteria->tmmcInit();
  feasst::WLTMMC mc(pair, criteria);
  mc.weight = 0.4;
  feasst::addTrialTransform(&mc,
    {{"transType", "translate"},
     {"maxMoveParam", "0.1"}});
  feasst::addTrialTransform(&mc,
    {{"transType", "rotate"},
     {"maxMoveParam", "0.1"}});

  ASSERT(nMolMax % 2 == 0, "assumes equimolar");
  mc.nMolSeek(nMolMax/2, ssA.str().c_str());
  mc.nMolSeek(nMolMax, ssB.str().c_str());

  mc.weight = 0.001;
  feasst::betaTrial(&mc, criteria->mBin());

  // geometric cluster algorithm
  mc.weight = 1./double(nMolMax);
  feasst::gcaTrial(&mc);

  // rigid cluster moves
  space->addTypeForCluster(0);
  space->addTypeForCluster(3);

  // Set the size of "premicellar aggregates"
  // This will print the concentration of "clusters" which are of size <= this value
  space->preMicellarAgg(5);

  mc.weight = 1./5./double(nMolMax);
  feasst::clusterTrial(&mc, "clustertrans");
  feasst::clusterTrial(&mc, "clusterrotate");

  // Note that the movie file outputs "patches" as inscribed spheres.
  // From PairPatchKFMulti::printxyz line 179:
  // make a patch by inscribing sphere of radius r2 inside bead, radius r1=sig
  //  distance between center of sphere and bead is a
  //  for given patch angle, use law of cosines to derive r2
  //  constraint r1 + eps = r2 + a, where eps is small, to avoid clipping
  const int nPrint = 1e4;
  mc.initLog("log", nPrint);
  mc.initColMat("colMat", nPrint);
  mc.initMovie("movie", nPrint);
  mc.initRestart("tmp/rst", nPrint);
  mc.setNFreqCheckE(nPrint, 1e-6);
  mc.setNFreqTune(nPrint);

  // reset cluster statistics after equilibration
  space->clusterReset();

  mc.runNumSweeps(10, 1e6);
}
