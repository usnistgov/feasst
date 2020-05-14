/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_fermi_jagla.h"

namespace feasst {

PairFermiJagla::PairFermiJagla(Space* space, const argtype &args)
  : Pair(space, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);
  argparse_.checkAllArgsUsed();
}

PairFermiJagla::PairFermiJagla(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
}

void PairFermiJagla::defaultConstruction_() {
  className_.assign("PairFermiJagla");
}

void PairFermiJagla::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double r2 = dx*dx + dy*dy + dz*dz;
  const double epsij = epsij_[iSiteType][jSiteType];
  const double sigij = sigij_[iSiteType][jSiteType];
  const double rinvsig = sqrt(r2)/sigij;
  *neighbor = 1;
  if (rinvsig <= r_s) {
    *energy = NUM_INF;
  } else {
    const double en_a = A_0/(1. + exp(A_1*(rinvsig - A_2)));
    const double en_b = -epsij/(1. + exp(B_1*(rinvsig - B_2)));
    const double en_c = epsilon_c*pow(sigma_c/(rinvsig - r_s), n_c);
    *energy = en_a + en_b + en_c;
  }
}

PairFermiJagla* PairFermiJagla::clone(Space* space) const {
  PairFermiJagla* p = new PairFermiJagla(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairFermiJagla> makePairFermiJagla(Space* space,
  const argtype &args) {
  return make_shared<PairFermiJagla>(space, args);
}

}  // namespace feasst

