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

int main() {  // SquareWell, SEMIGRAND_EXPANDED
  feasst::ranInitByDate();
  auto space = feasst::makeSpace(
   {{"dimen", "3"},
    {"boxLength", "8."}});
  space->equiMolar(1);
  feasst::PairSquareWell pair(space.get(), {{"rCut", "1.5"}});
  std::stringstream addMolTypeA, addMolTypeB;
  addMolTypeA << space->install_dir() << "/forcefield/data.lj";
  addMolTypeB << space->install_dir() << "/forcefield/data.ljb";
  pair.initData(addMolTypeA.str().c_str());
  pair.initData(addMolTypeB.str().c_str());
  pair.epsijset(0, 0, 1.);
  pair.epsijset(0, 1, sqrt(0.93));
  pair.epsijset(1, 1, 0.93);
  pair.rCutijset(0, 0, 1.5);
  pair.rCutijset(0, 1, 1.5);
  pair.rCutijset(1, 1, 1.5);
  pair.initEnergy();

  // initialize acceptance criteria
  auto criteria = feasst::makeCriteriaWLTMMC(
   {{"beta", feasst::str(1./5.)},
    {"activ", feasst::str(exp(-2))},  // activity of A
    {"mType", "nmol"},
    {"nMax", "30"}
   });
  criteria->addActivity(exp(-2));  // activity of B
  criteria->collectInit();
  criteria->tmmcInit();

  // initialize MC simulation object
  feasst::WLTMMC mc(space.get(), &pair, criteria.get());
  mc.weight = 1.;
  feasst::transformTrial(&mc, "translate");
  mc.weight = 1/4./4.;
  feasst::insertDeleteTrial(&mc, addMolTypeA.str().c_str());
  feasst::insertDeleteTrial(&mc, addMolTypeB.str().c_str());
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
//  mc.initWindows(2.0, 0);
  mc.runNumSweeps(1, -1);
}


