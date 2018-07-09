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

int main() {  // SPCE, SRSW_NVTMC_SiteCut
  feasst::ranInitByDate();
  auto space = feasst::makeSpace(
    {{"dimen", "3"},
     {"boxLength", "24.8586887"}});

  auto pair = feasst::makePairLJCoulEwald(space,
    {{"rCut", feasst::str(space->minl()/2.)},
     {"molTypeInForcefield", "data.spce"},
     {"alphaL", "5.6"},
     {"k2max", "38"}});
  pair->initAtomCut(1);
  pair->equateRcutForAllTypes();

  // acceptance criteria
  const double temperature = 298;  // Kelvin
  auto criteria = feasst::makeCriteriaMetropolis(
    {{"beta", feasst::str(1./(temperature*feasst::idealGasConstant/1e3))},  // mol/KJ
     {"activ", "1."}});
  feasst::MC mc(pair, criteria);
  feasst::transformTrial(&mc, "translate", 0.1);
  feasst::transformTrial(&mc, "rotate", 0.1);
  const int nPrint = 1e4;
  mc.initLog("log", nPrint);
  mc.initMovie("movie", nPrint);
  mc.initRestart("tmp/rst", nPrint);
  mc.setNFreqTune(nPrint);
  mc.setNFreqCheckE(nPrint, 1e-6);
  mc.nMolSeek(512);
  mc.runNumTrials(1e6);   // run equilibration

  // Run the production simulation
  mc.initProduction();
  mc.zeroStat();
  mc.setNFreqTune(0);
  mc.runNumTrials(1e6);

  // Check average energy against Gerhard Hummer
  // https://doi.org/10.1063/1.476834
  const double
    peAv = mc.peAccumulator().average()/static_cast<double>(space->nMol()),
    peStd = mc.peAccumulator().blockStdev()/static_cast<double>(space->nMol()),
    pePublish = -46.82,  // published value
    pePublishStd = 0.02;
  ASSERT(fabs(peAv - pePublish) < 2.576*(pePublishStd+peStd),  // 99% confidence
    "ERROR: The average potential energy(" << peAv << " +/- " << peStd
    << ") did not match the published value (" << pePublish << " +/- " << pePublishStd << ")");
}
