/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#include "./pair_squarewell.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairSquareWell::PairSquareWell(Space* space, const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairSquareWell::PairSquareWell(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
}

void PairSquareWell::defaultConstruction_() {
  className_.assign("PairSquareWell");
}

void PairSquareWell::initHardSphere(const int itype, const int jtype) {
  rCutijset(itype, jtype, sigij_[itype][jtype]);
}

void PairSquareWell::multiPartEnerAtomCutInner(const double &r2,
  const int &itype,
  const int &jtype) {
  const double sigSq = sigij_[itype][jtype]*sigij_[itype][jtype];
  if (r2 < sigSq) {
    // hard sphere if less than sigma
    peSRone_ += NUM_INF;
  } else {
    // otherwise, square well
    peSRone_ -= epsij_[itype][jtype];
  }
}

PairSquareWell* PairSquareWell::clone(Space* space) const {
  PairSquareWell* p = new PairSquareWell(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairSquareWell> makePairSquareWell(Space* space, const double rCut) {
  return make_shared<PairSquareWell>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

