/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "pair_hard_sphere.h"
#include "mc.h"
#include "trial_transform.h"

int main() {  // HardSphere, NVTMC
  feasst::Space space(3);
  const double rho = 1e-3;  // number density
  const int nMol = 500;     // number of particles
  space.lset(pow(double(nMol)/rho, 1./3.));   // set the cubic PBCs
  feasst::PairHardSphere pair(&space);
  stringstream molNameSS;
  molNameSS << space.install_dir() << "/forcefield/data.atom";
  pair.initData(molNameSS.str().c_str());
  pair.initEnergy();
  feasst::CriteriaMetropolis criteria(1., 1.);
  feasst::MC mc(&space, &pair, &criteria);
  feasst::transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(nMol);
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.initRestart("tmp/rst", 1e4);
  mc.setNFreqTune(1e4);
  mc.runNumTrials(1e7);   // run equilibration

  // Run the production simulation
  mc.runNumTrials(1e7);
}
