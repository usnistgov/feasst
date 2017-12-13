/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_hard_sphere.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairHardSphere::PairHardSphere(Space* space)
  : Pair(space) {
  defaultConstruction_();
}

PairHardSphere::PairHardSphere(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
}

void PairHardSphere::defaultConstruction_() {
  className_.assign("PairHardSphere");
}

void PairHardSphere::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double sigij = sigij_[iSiteType][jSiteType];
  const double r2 = dx*dx + dy*dy + dz*dz;
  if (r2 < sigij*sigij) {
    *energy = NUM_INF;
    *neighbor = 1;
    return;
  }
  *energy = 0;
  *neighbor = 0;
}

void PairHardSphere::initPairParam(const vector<double> eps,
  const vector<double> sig, const vector<double> sigref) {
  Pair::initPairParam(eps, sig, sigref);
  sig2rCut();
  rCut_ = rCutMaxAll_;
}

PairHardSphere* PairHardSphere::clone(Space* space) const {
  PairHardSphere* p = new PairHardSphere(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairHardSphere> makePairHardSphere(Space* space) {
  return make_shared<PairHardSphere>(space);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

