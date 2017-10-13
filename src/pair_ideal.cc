/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./pair_ideal.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairIdeal::PairIdeal(Space* space, const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairIdeal::PairIdeal(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
}

void PairIdeal::defaultConstruction_() {
  className_.assign("PairIdeal");
}

void PairIdeal::multiPartEnerAtomCutInner(const double &r2,
  const int &itype,
  const int &jtype) {
  if (static_cast<int>(r2) < itype*jtype) {}  // remove unused variable warning
}

PairIdeal* PairIdeal::clone(Space* space) const {
  PairIdeal* p = new PairIdeal(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairIdeal> makePairIdeal(Space* space, const double rCut) {
  return make_shared<PairIdeal>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

