/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "pair.h"
#include "mc.h"
#include "trial_transform.h"

using namespace feasst;

// Define a new Pair class to simulate the Jagla potential
// sigma is the hard particle size and rCut is the extent of the subsequent
// linear ramp. Epsilon is the size of the repulsion at r=sigma.
class PairJagla : public Pair {
 public:
  PairJagla(Space * space, const double rCut) : Pair(space, rCut) {}
  ~PairJagla() {}

  // Overloaded virtual function from pair.h
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype) {
    const double sigij = sigij_[itype][jtype];
    if (r2 < sigij*sigij) {
      peSRone_ += NUM_INF;
    } else {
      const double z = (sqrt(r2) - sigij) / (rCutij_[itype][jtype] - sigij);
      peSRone_ += epsij_[itype][jtype]*z;
    }
  }
};

int main() {  // JAGLA, NVTMC
  Space space(3);
  const double rho = 1e-2;  // number density
  const int nMol = 500;     // number of particles
  space.lset(pow(double(nMol)/rho, 1./3.));   // set the cubic PBCs
  stringstream molNameSS;
  molNameSS << space.install_dir() << "/forcefield/data.lj";
  space.addMolInit(molNameSS.str().c_str());
  PairJagla pair(&space, 2);   // potential truncation at 3
  pair.initData(molNameSS.str());
  pair.rCutijset(0, 0, pair.rCut());
  const double temperature = 1.;
  CriteriaMetropolis criteria(1./temperature, 1.);
  MC mc(&space, &pair, &criteria);
  transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(nMol, molNameSS.str().c_str());
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  mc.initRestart("tmp/rst", 1e4);
  mc.setNFreqTune(1e4);
  mc.runNumTrials(1e7);   // run equilibration

  // Run the production simulation
  mc.runNumTrials(1e7);
}
