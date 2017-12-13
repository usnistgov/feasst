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

int main() {  // LJ, NPT
  feasst::ranInitByDate();
  feasst::Space space(3);
  space.initBoxLength(8);
  feasst::PairLJ pair(&space, {{"rCut", "3."}});
  std::stringstream ss;
  ss << space.install_dir() << "/forcefield/data.lj";
  pair.initData(ss.str());
  pair.cutShift(1);
  pair.initEnergy();
  feasst::CriteriaMetropolis crit(1./1.5, exp(-4));
  feasst::MC mc(&space, &pair, &crit);
  feasst::transformTrial(&mc, "translate", 1);
  mc.nMolSeek(50);
  crit.pressureset(0.1);
  feasst::transformTrial(&mc, "lxmod", 0.001);
  feasst::transformTrial(&mc, "lymod", 0.001);
  feasst::transformTrial(&mc, "lzmod", 0.001);
  // feasst::transformTrial(&mc, "vol", 0.001);
  const int nfreq = 1e3, ncfreq = 1e4;
  mc.initLog("log", nfreq);
  mc.setNFreqCheckE(ncfreq, 2e-1);
  mc.setNFreqTune(ncfreq);
  mc.initMovie("movie", ncfreq);
  mc.initRestart("tmp/rst", ncfreq);
  mc.runNumTrials(20*ncfreq);
}
