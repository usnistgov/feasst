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

namespace feasst {

PairHardSphere::PairHardSphere(Space* space)
  : Pair(space) {
  defaultConstruction_();
}

PairHardSphere::PairHardSphere(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  string strtmp = fstos("reduction_factor", fileName);
  if (!strtmp.empty()) {
    setReductionFactor(stod(strtmp));
  }
}

void PairHardSphere::defaultConstruction_() {
  className_.assign("PairHardSphere");
  setReductionFactor();
}

void PairHardSphere::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double sigij = reduction_factor_*sigij_[iSiteType][jSiteType];
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

void PairHardSphere::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  if (reduction_factor_ != 1.) {
    file << "# reduction_factor " << reduction_factor_ << endl;
  }
}

shared_ptr<PairHardSphere> makePairHardSphere(Space* space) {
  return make_shared<PairHardSphere>(space);
}

}  // namespace feasst

