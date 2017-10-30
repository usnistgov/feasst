/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pagespace.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "pair_squarewell.h"
#include "mc_wltmmc.h"
#include "trial_transform.h"
#include "ui_abbreviated.h"
#include "trial_swap.h"
#include "trial_xswap.h"

int main() {  // SquareWell, SEMIGRAND
  feasst::ranInitByDate();
  feasst::Space space(3, 0);
  space.lset(8.);
  feasst::PairSquareWell pair(&space, 2.);
  std::stringstream addMolTypeA, addMolTypeB;
  addMolTypeA << space.install_dir() << "/forcefield/data.lj";
  addMolTypeB << space.install_dir() << "/forcefield/data.ljb";
  pair.initData(addMolTypeA.str().c_str());
  pair.initData(addMolTypeB.str().c_str());
  pair.rCutijset(0, 0, 2.);
  pair.rCutijset(0, 1, 1.75);
  pair.rCutijset(1, 1, 1.5);
  pair.initEnergy();

  // initialize acceptance criteria
  const double beta = 1./5., activ = exp(-4.568214);  // activity of A
  const int nMolMax = 30, nMolMin = 0;
  feasst::CriteriaWLTMMC crit(beta, activ, "nmol", nMolMin, nMolMax);
  crit.addActivity(exp(-3));  //!< activity of B

  // initialize MC simulation object
  feasst::WLTMMC mc(&space, &pair, &crit);
  mc.weight=3./4.;
  feasst::transformTrial(&mc, "translate");
  mc.weight=1./16.;
  feasst::insertDeleteTrial(&mc, addMolTypeA.str().c_str());
  feasst::insertDeleteTrial(&mc, addMolTypeB.str().c_str());
  feasst::swapTrial(&mc, addMolTypeA.str().c_str(), addMolTypeB.str().c_str());
  feasst::xswapTrial(&mc);

  // output log, lnpi and movie
  const int nfreq = 1e4, ncfreq = 1e6;
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(nfreq);
  mc.initMovie("movie", nfreq);
  //mc.initXTC("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  //production tmmc simulation
  crit.collectInit();
  crit.tmmcInit();
  mc.runNumSweeps(1, -1);
}


