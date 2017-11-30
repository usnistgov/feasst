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

int main() {  // LJ, MULTISITE
  feasst::ranInitByDate();
  feasst::Space space(3);
  space.initBoxLength(9.);
  stringstream addMolType;
  addMolType << space.install_dir() << "/forcefield/data.cg3_60_1_1";
  feasst::PairLJ pair(&space, 3.,  // potential cut-off
    {{"molType", addMolType.str()}});
  pair.rCutijset(1, 1, pair.rCut());
  pair.linearShiftijset(1, 1, 1);
  pair.initWCA(1, 2);
  pair.initWCA(2, 2);
  pair.initEnergy();

  // Initialize Monte Carlo acceptance and moves
  const double beta = 1./0.225;
  const double activ = exp(-6.);
  const int nMolMin = 0, nMolMax = 30;
  feasst::CriteriaWLTMMC crit(beta, activ, "nmol", nMolMin - 0.5, nMolMax + 0.5, nMolMax - nMolMin + 1);
  feasst::WLTMMC mc(&space, &pair, &crit);
  mc.weight = 1;
  feasst::transformTrial(&mc, "translate");
  feasst::transformTrial(&mc, "rotate");
  mc.weight = 1./4.;
  feasst::insertDeleteTrial(&mc, addMolType.str().c_str());

  // print to files
  const int nfreq = 1e4, ncfreq = 1e5;
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(ncfreq);
  mc.initMovie("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  // Turn on Wang-Landau and transition-matrix
  crit.collectInit();
  crit.tmmcInit();

  mc.runNumSweeps(1);
}
