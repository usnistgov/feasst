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

// Define a new Analyze class in order to compute the average potential energy
class AnalyzeMonkeyPatch : public feasst::Analyze {
 public:
  AnalyzeMonkeyPatch(feasst::Pair *pair) : Analyze(pair) {}
  ~AnalyzeMonkeyPatch() {}

  // Access the average and block std from the accumulator
  feasst::Accumulator pe;

  // Overwrite the virtual function to accumulate intensive potential energy
  void update() {
    pe.accumulate(pair_->peTot()/double(space()->nMol()));
  }
  void write() {
    cout << pe.average() << " +/- " << pe.blockStdev() << endl;
  }
};

int main() {  // LJ, SRSW_NVTMC
  feasst::Space space(3);
  const double rho = 1e-3;  // number density
  const int nMol = 500;     // number of particles
  space.lset(pow(double(nMol)/rho, 1./3.));   // set the cubic PBCs
  stringstream molNameSS;
  molNameSS << space.install_dir() << "/forcefield/data.lj";
  space.addMolInit(molNameSS.str().c_str());
  feasst::PairLJ pair(&space, 3);   // potential truncation at 3
  pair.initEnergy();
  const double temperature = 0.9;
  feasst::CriteriaMetropolis criteria(1./temperature, 1.);
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(nMol);
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.initRestart("tmp/rst", 1e4);
  mc.setNFreqTune(1e4);
  mc.runNumTrials(1e7);   // run equilibration

  // Initialize the custom analysis to compute average energy
  shared_ptr<AnalyzeMonkeyPatch> an =
    make_shared<AnalyzeMonkeyPatch>(&pair);
  an->initFreq(1);    // frequency that Analyze::update() is called
  an->initPrintFreq(1e7);
  mc.initAnalyze(an);

  // Run the production simulation
  mc.runNumTrials(1e7);

  // Check average energy against the NIST SRSW
  // https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
  // https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website
  const double
    peAv = an->pe.average(),
    peStd = an->pe.blockStdev(),
    peSRSW = -9.9165E-03,
    peSRSWstd = 1.89E-05;
  ASSERT(fabs(peAv - peSRSW) < 2.576*(peSRSWstd+peStd),  // 99% confidence
    "ERROR: The average potential energy(" << peAv << " +/- " << peStd
    << ") did not match the SRSW (" << peSRSW << " +/- " << peSRSWstd << ")");
}
