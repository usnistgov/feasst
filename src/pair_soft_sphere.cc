/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_soft_sphere.h"

namespace feasst {

PairSoftSphere::PairSoftSphere(Space* space, const argtype &args)
  : Pair(space, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);
  argparse_.checkAllArgsUsed();
}

PairSoftSphere::PairSoftSphere(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  string str = fstos("npotentialparameter", fileName);
  if (!str.empty()) {
    initExponent(stoi(str));
  }
}

void PairSoftSphere::defaultConstruction_() {
  className_.assign("PairSoftSphere");
  initExponent();
}

void PairSoftSphere::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# npotentialparameter " << n_ << endl;
}

void PairSoftSphere::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double rCut = rCutij_[iSiteType][jSiteType];
  const double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 < rCut*rCut) {
    const double sigij = sigij_[iSiteType][jSiteType];
    *energy = pow(sigij*sigij/r2, n_/2.);
    *neighbor = 1;
    return;
  }
  *energy = 0;
  *neighbor = 0;
  return;
}

double PairSoftSphere::b2reduced() {
  ASSERT(n_ == 12, "virial coefficient hard coded for n=12 for gamma[1-3/n]");
  const double gamma = 1.225416702465178;
  const double sigij = sigij_[0][0];
  return 2.*PI/3.*pow(sigij, 3)*gamma;
}

PairSoftSphere* PairSoftSphere::clone(Space* space) const {
  PairSoftSphere* p = new PairSoftSphere(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairSoftSphere> makePairSoftSphere(Space* space,
  const argtype &args) {
  return make_shared<PairSoftSphere>(space, args);
}

}  // namespace feasst

