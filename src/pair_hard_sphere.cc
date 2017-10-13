/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./pair_hard_sphere.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairHardSphere::PairHardSphere(Space* space, const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairHardSphere::PairHardSphere(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
}

void PairHardSphere::defaultConstruction_() {
  className_.assign("PairHardSphere");
}

void PairHardSphere::multiPartEnerAtomCutInner(const double &r2,
  const int &itype,
  const int &jtype) {
  const double sigij = sigij_[itype][jtype];
  if (r2 < sigij*sigij) {
    peSRone_ += NUM_INF;
  }
  // cout << "in Inner: r2 " << r2 << " itype " << itype << " jtype " << jtype
  //   << " peSRone " << peSRone_ << " sigij " << sigij << endl;
}

PairHardSphere* PairHardSphere::clone(Space* space) const {
  PairHardSphere* p = new PairHardSphere(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairHardSphere> makePairHardSphere(Space* space, const double rCut) {
  return make_shared<PairHardSphere>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

