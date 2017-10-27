/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_hard_circle.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

PairHardCircle::PairHardCircle(Space* space, const double rCut)
  : Pair(space, rCut) {
  defaultConstruction_();
}

PairHardCircle::PairHardCircle(Space* space, const char* fileName)
  : Pair(space, fileName) {
  defaultConstruction_();
  dCircle_ = fstod("dCircle", fileName);
  rDep_ = fstod("rDep", fileName);
}

void PairHardCircle::defaultConstruction_() {
  className_.assign("PairHardCircle");
  dCircle_ = 1;
  initRDep();
  ASSERT(space_->dimen() == 2,
    "Hard circle potential requires spatial dimenion(" << space_->dimen()
    << ") of 2");
}

void PairHardCircle::writeRestart(const char* fileName) {
  writeRestartBase(fileName);
  std::ofstream file(fileName, std::ios_base::app);
  file << "# dCircle " << dCircle_ << endl;
  file << "# rDep " << rDep_ << endl;
}

void PairHardCircle::multiPartEnerAtomCutInner(const double &r2,
  const int &itype,
  const int &jtype) {
  if (itype == jtype) {}  // remove unused parameter warning
  // hard sphere
  if (r2 < dCircle_*dCircle_) {
    peSRone_ += NUM_INF;
  } else {
    const double r = sqrt(r2),
                 R = 0.5*dCircle_ + rDep_;
    peSRone_ -= (2*R*R*acos(r*0.5/R) - r*sqrt(R*R-r*r*0.25))/PI/rDep_/rDep_;
  }
}

PairHardCircle* PairHardCircle::clone(Space* space) const {
  PairHardCircle* p = new PairHardCircle(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairHardCircle> makePairHardCircle(Space* space, const double rCut) {
  return make_shared<PairHardCircle>(space, rCut);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

