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

int main() {  // SquareWell, FLOPPY
  feasst::ranInitByDate();
  feasst::Space space(3);
  space.initBoxLength(8.);
  space.setXYTilt(1.5);
  space.setXZTilt(1.5);
  space.setYZTilt(1.5);
  feasst::PairSquareWell pair(&space, 2.);
  std::stringstream addMolType;
  addMolType << space.install_dir() << "/forcefield/data.lj";
  pair.initData(addMolType.str());
  pair.rCutijset(0, 0, pair.rCut());
  pair.initEnergy();
  feasst::CriteriaMetropolis crit(0.011, exp(-4));
  feasst::MC mc(&space, &pair, &crit);
  mc.nMolSeek(100);  // important to set this before pressure set
  feasst::transformTrial(&mc, "translate", 1);
  feasst::transformTrial(&mc, "xytilt", 0.01);
  feasst::transformTrial(&mc, "xztilt", 0.01);
  feasst::transformTrial(&mc, "yztilt", 0.01);
  crit.pressureset(100);
  feasst::transformTrial(&mc, "lxmod", 0.1);
  feasst::transformTrial(&mc, "lymod", 0.1);
  feasst::transformTrial(&mc, "lzmod", 0.1);
  const int nfreq = 1e3, ncfreq = 1e4;
  mc.initLog("log", nfreq);
  mc.setNFreqCheckE(ncfreq, 2e-1);
  mc.initMovie("movie",ncfreq);
  mc.initRestart("tmp/rst",    ncfreq);
  mc.runNumTrials(20*ncfreq);
}
