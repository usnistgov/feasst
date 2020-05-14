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

// Define a new Analyze class in order to compute the average potential energy
class AnalyzePE_TC13 : public feasst::Analyze {
 public:
  AnalyzePE_TC13(shared_ptr<feasst::Pair> pair) : Analyze(pair) {}
  feasst::Accumulator pe;  // average and block std from the accumulator
  void update() {          // overwrite the virtual function to update pe
    pe.accumulate(pair_->peTot()/double(space()->nMol()));
  }
  void write() {           // overwrite to print pe
    cout << "Potential energy: "
         << pe.average() << " +/- " << pe.blockStdev() << endl;
  }
};

int main() {  // LJ, SRSW_NVTMC
  feasst::ranInitByDate();
  const double rho = 1e-3;  // number density
  const int nMol = 500;     // number of particles
  auto space = feasst::makeSpace(
    {{"dimen", "3"},
     {"boxLength", feasst::str(pow(double(nMol)/rho, 1./3.))}});
  auto pair = feasst::makePairLJ(space,
    {{"rCut", "3"},
     {"cutType", "lrc"},
     {"molTypeInForcefield", "data.lj"}});
  auto criteria = feasst::makeCriteriaMetropolis(
    {{"beta", feasst::str(1./0.9)}});   // beta=1/kT
  feasst::MC mc(pair, criteria);
  feasst::addTrialTransform(&mc,
    {{"transType", "translate"},
     {"maxMoveParam", "0.1"}});
  mc.nMolSeek(nMol);
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.initRestart("tmp/rst", 1e4);
//  mc.setNFreqTune(1e4);
  mc.setNFreqCheckE(1e4, 1e-6);

//  auto xtc = feasst::makeAnalyzeTRAJ(pair.get(),
//    {{"nFreqPrint", feasst::str(1e4)},
//     {"fileName", "movie2"},
//     {"format", "xtc"}});
//  mc.initAnalyze(xtc);

  mc.runNumTrials(1e5);   // run equilibration

  // Initialize the custom analysis to compute average energy
  auto an = make_shared<AnalyzePE_TC13>(pair);
  an->initFreq(1);    // frequency that Analyze::update() is called
  an->initFreqPrint(1e5);
  mc.initProduction();
  mc.initAnalyze(an);

  //mc.setWallClockMaxInHours(1./60./60.);

  // Run the production simulation
  mc.runNumTrials(1e5);

  // Check average energy against the NIST SRSW
  // https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm
  // https://www.nist.gov/programs-projects/nist-standard-reference-simulation-website
}
