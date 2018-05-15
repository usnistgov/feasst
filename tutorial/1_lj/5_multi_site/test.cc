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
  auto space = feasst::makeSpace(
      {{"dimen", "3"},
       {"boxLength", "9"}});
  auto pair = feasst::makePairLJ(space,
    {{"rCut", "3."},
     {"molTypeInForcefield", "data.cg3_60_1_1"}});
  pair->rCutijset(1, 1, pair->rCut());
  pair->linearShiftijset(1, 1, 1);
  pair->initWCA(1, 2);
  pair->initWCA(2, 2);
  pair->initEnergy();

  // Use Wang-Landau and transition-matrix
  auto crit = feasst::makeCriteriaWLTMMC(
    {{"beta", feasst::str(1./0.225)},
     {"activ", feasst::str(exp(-6))},
     {"mType", "nmol"},
     {"nMax", "30"}});
  crit->collectInit();
  crit->tmmcInit();

  // Initialize Monte Carlo acceptance and moves
  feasst::WLTMMC mc(pair, crit);
  mc.weight = 1;
  feasst::transformTrial(&mc, "translate");
  feasst::transformTrial(&mc, "rotate");
  mc.weight = 1./4.;
  feasst::insertDeleteTrial(&mc, space->addMolListType(0).c_str());

  // print to files
  const int nfreq = 1e4, ncfreq = 1e5;
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(ncfreq);
  mc.initMovie("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  mc.runNumSweeps(1);
}
