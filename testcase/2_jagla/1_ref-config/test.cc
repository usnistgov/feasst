/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
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

int main() {  // JAGLA, REF_CONFIG
  Space space(3);
  stringstream molNameSS;
  PairJagla pair(&space, 2);   // potential truncation represents end of "ramp"
  molNameSS << space.install_dir() << "/forcefield/data.atom";

  // read from the data file in order to set the hard particle size (sigma)
  // and the repulsive potential value at sigma (epsilon)
  pair.initData(molNameSS.str());

  // set the cut-off to the end of the "ramp" potential
  pair.rCutijset(0, 0, pair.rCut());

  // first, test the interaction between two particles
  vector<double> xAdd(space.dimen(), 0.);  // position to add is origin
  space.xAdd = xAdd;        // tell space that the next addMol goes to xAdd
  pair.addMol();            // add first molecule (randomly if xAdd not set)
  const double r = 1.2345;
  xAdd[0] = r;              // position of second particle
  space.xAdd = xAdd;        // xAdd must be reset after every addMol
  pair.addMol();            // add second particle
  pair.initEnergy();        // compute energy
  const double peExpected = r - 1.;
  ASSERT(fabs(pair.peTot() - peExpected) < feasst::DTOL,
    "analytical comparison");
}
