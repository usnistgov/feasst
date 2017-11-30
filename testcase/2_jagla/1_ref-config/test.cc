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

// Define a new Pair class to simulate the Jagla potential
// sigma is the hard particle size and rCut is the extent of the subsequent
// linear ramp. Epsilon is the size of the repulsion at r=sigma.
// Non-dimensionalize by sigma and epsilon such that the only parameter
// is now k_B T / epsilon and Delta = (rCut-sigma)/sigma = rCut/sigma - 1
class PairJagla : public feasst::Pair {
 public:
  PairJagla(feasst::Space * space, const double rCut)
    : feasst::Pair(space, rCut) {}
  ~PairJagla() {}

  // Overloaded virtual function from pair.h
  void pairSiteSite_(
    const int &iSiteType,  //!< type of first site
    const int &jSiteType,  //!< type of second site
    double * energy,      //!< energy of interaction
    double * force,       //!< force/rij of interaction
    int * neighbor,       //!< 1 if neighbor, 0 otherwise
    const double &dx,      //!< x-dimension separation
    const double &dy,      //!< y-dimension separation
    const double &dz) {    //!< z-dimension separation
    *neighbor = 1;
    const double r2 = dx*dx + dy*dy + dz*dz;
    // hard sphere if less than sigma
    if (r2 < 1.) {
      *energy = NUM_INF;
    // otherwise, linear ramp
    } else {
      *energy += (sqrt(r2) - 1.)/(rCutij_[iSiteType][jSiteType] - 1.);
    }
  }
};

int main() {  // JAGLA, REF_CONFIG
  feasst::Space space(3);
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
