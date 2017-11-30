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

int main() {  // LJ, SRSW_REFCONF
  feasst::Space space(3);
  space.initBoxLength(8);
  feasst::PairLJ pair(&space, 3.,  // potential truncation
    {{ "cutType", "lrc" }});            // use long range corrections

  // create clones of Space and Pair to perform two separate tests
  feasst::Space * space2 = space.clone();
  feasst::PairLJ * pair2 = pair.clone(space2);

  // first, test the interaction between two particles
  vector<double> xAdd(space.dimen(), 0.);   // position to add is origin
  pair.addMol(xAdd);                        // add first molecule at xAdd
  const double r = 1.2345;
  xAdd[0] = r;              // position of second particle
  pair.addMol(xAdd);        // add second particle
  pair.initEnergy();        // compute energy
  ASSERT(fabs(pair.peLJ() - 4*(pow(r, -12) - pow(r, -6))) < feasst::DTOL,
    "LJ equation failed");
  ASSERT(fabs(pair.peLRC() - (8./3.)*feasst::PI*pow(space.nMol(), 2)/space.vol()
    *((1./3.)*pow(pair.rCut(), -9) - pow(pair.rCut(), -3))) < feasst::DTOL,
    "LRC equation failed");

  // second, compare with the reference configuration
  stringstream config;
  config << space.install_dir() <<
    "/testcase/1_lj/1_ref-config/lj_sample_config_periodic4.xyz";
  std::ifstream file(config.str().c_str());   // load the xyz file
  pair2->readxyz(file);                       // read the xyz file
  pair2->initEnergy();                        // compute energy
  const double peLJ = -16.790321304625856, peLRC = -0.5451660014945704;
  ASSERT(fabs(pair2->peLRC() - peLRC) < feasst::DTOL, "SRSW LJ-LRC comparison failed");
  ASSERT(fabs(pair2->peLJ() - peLJ) < feasst::DTOL, "SRSW LJ comparison failed");
}
