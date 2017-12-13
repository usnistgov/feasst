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

PairHardCircle::PairHardCircle(Space* space, const argtype &args)
  : Pair(space, args) {
  defaultConstruction_();
  argparse_.initArgs(className_, args);
  argparse_.checkAllArgsUsed();
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

void PairHardCircle::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double r2 = dx*dx + dy*dy + dz*dz;
  *energy = 0;
  *neighbor = 0;
  if (r2 < dCircle_*dCircle_) {
    // hard sphere
    *energy = NUM_INF;
    *neighbor = 1;
  } else {
    const double r = sqrt(r2),
                 R = 0.5*dCircle_ + rDep_;
    *energy = -(2*R*R*acos(r*0.5/R) - r*sqrt(R*R-r*r*0.25))/PI/rDep_/rDep_;
    *neighbor = 1;
  }
}

PairHardCircle* PairHardCircle::clone(Space* space) const {
  PairHardCircle* p = new PairHardCircle(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairHardCircle> makePairHardCircle(Space* space,
  const argtype &args) {
  return make_shared<PairHardCircle>(space, args);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

