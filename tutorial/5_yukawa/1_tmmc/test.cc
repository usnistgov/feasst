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

int main() {  // Yukawa, NVTMC
  feasst::ranInitByDate();
  feasst::Space space(3);
  space.initBoxLength(12.);
  //space.updateCells(rCut);
  std::stringstream addMolType;
  addMolType << space.install_dir() << "/forcefield/data.atom";
  feasst::PairLJ pair(&space, {{"rCut", "5"}, {"molType", addMolType.str()}});
  pair.initExpType(3);

  // Set the Yukawa potential
  // A*exp(-K r)/r; in bolhuis, yukawa = A zetaexp(-r/zeta)/r... Thus, A=A zeta, K = 1/zeta
  pair.initScreenedElectro(2*1.794, 1./1.794);
  pair.cutShift(1);
  pair.initEnergy();

  const int nMolMin = 0, nMolMax = 30;
  feasst::CriteriaWLTMMC crit(7., exp(-2), "nmol", nMolMin, nMolMax);
  feasst::WLTMMC mc(&space, &pair, &crit);
  feasst::transformTrial(&mc, "translate");
  mc.weight = 1./4.;
  feasst::deleteTrial(&mc);
  feasst::addTrial(&mc, addMolType.str().c_str());

  // output log, lnpi and movie
  const int nfreq = 1e5, ncfreq = 1e6;
  mc.initLog("log", nfreq);
  mc.initColMat("colMat", ncfreq);
  mc.setNFreqCheckE(ncfreq, 2e-4);
  mc.setNFreqTune(ncfreq);
  mc.initMovie("movie", nfreq);
  mc.initXTC("movie", nfreq);
  mc.initRestart("tmp/rst", ncfreq);

  //production tmmc simulation
  mc.nMolSeek(nMolMin, 1e6);
  crit.collectInit();
  crit.tmmcInit();
  mc.runNumSweeps(1, -1);
}


