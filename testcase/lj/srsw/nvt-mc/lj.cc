/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "pair_lj.h"
#include "mc.h"
#include "trial_transform.h"

using namespace feasst;

class AnalyzeMonkeyPatch : public Analyze {
 public:
  AnalyzeMonkeyPatch(Pair *pair) : Analyze(pair) {}
  ~AnalyzeMonkeyPatch() {}
  Accumulator pe;
  void update() {
    pe.accumulate(pair_->peTot()/double(space()->nMol()));
  }
  void write() {
    cout << pe.average() << " +/- " << pe.blockStdev() << endl;
  }
};

int main() {
  Space space(3, 0);
  const double rho = 1e-3;
  const int nMol = 500;
  for (int dim = 0; dim < space.dimen(); ++dim) {
    space.lset(pow(double(nMol)/rho, 1./3.), dim);
  }
  stringstream molNameSS;
  molNameSS << space.install_dir() << "/forcefield/data.lj";
  space.addMolInit(molNameSS.str().c_str());
  PairLJ pair(&space, 3);   // potential truncation at 3
  pair.initEnergy();
  const double temperature = 0.9;
  CriteriaMetropolis criteria(1./temperature, 1.);
  MC mc(&space, &pair, &criteria);
  transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(nMol, molNameSS.str().c_str());
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.initRestart("tmp/rst", 1e4);
  mc.setNFreqTune(1e4);
  mc.runNumTrials(1e6);   // run equilibration
  shared_ptr<AnalyzeMonkeyPatch> an =
    make_shared<AnalyzeMonkeyPatch>(&pair);
  an->initFreq(1);
  an->initPrintFreq(1e5);
  mc.initAnalyze(an);
  mc.runNumTrials(1e6);

  // Check average energy against the NIST SRSW
  // https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
  // https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website
  const double
    peAv = an->pe.average(),
    peStd = an->pe.blockStdev(),
    peSRSW = -9.9165E-03,
    peSRSWstd = 1.89E-05;
  ASSERT(fabs(peAv - peSRSW) < peSRSWstd+peStd,
    "ERROR: The average potential energy(" << peAv << " +/- " << peStd
    << ") did not match the SRSW (" << peSRSW << " +/- " << peSRSWstd << ")");
}
