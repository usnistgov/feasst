/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "./pair_hs_coul_ewald.h"

namespace feasst {

PairHSCoulEwald::PairHSCoulEwald(Space* space, const argtype &args)
  : PairLJCoulEwald(space, args) {
  defaultConstruction_();
}

PairHSCoulEwald::PairHSCoulEwald(Space* space, const char* fileName)
  : PairLJCoulEwald(space, fileName) {
  defaultConstruction_();
}

void PairHSCoulEwald::defaultConstruction_() {
  className_.assign("PairHSCoulEwald");
  lrcFlag = 0;
  chargeConversion_ = 1.;
}

void PairHSCoulEwald::pairSiteSite_(const int &iSiteType, const int &jSiteType,
  double * energy, double * force, int * neighbor, const double &dx,
  const double &dy, const double &dz) {
  const double r2 = dx*dx + dy*dy + dz*dz;
  const double sigij = sigij_[iSiteType][jSiteType];
  *energy = 0.;
  *force = 0.;
  *neighbor = 1;
  if (r2 < sigij*sigij) {
    *energy += NUM_INF;
    peLJone_ += NUM_INF;
    *force += NUM_INF;
    peSRone_ += *energy;
    return;
  }

  // charge interactions
  const double enq = q_[iSiteType]*q_[jSiteType]*erft_.eval(r2);
  *energy += enq;
  peQRealone_ += enq;
  peSRone_ += enq;
  *force += q_[iSiteType]*q_[jSiteType]
            *(2.*alpha*exp(-alpha*alpha*r2)/sqrt(PI));
}

int PairHSCoulEwald::multiPartEnerAtomCutInner(const double &r2,
  const int &itype, const int &jtype) {
  const double sigij = sigij_[itype][jtype];
  double en = 0.;
  if (r2 < sigij*sigij) {
    en = NUM_INF;
  }
  peLJone_ += en;

  // charge interactions
  const double enq = q_[itype]*q_[jtype]*erft_.eval(r2);
  peQRealone_ += enq;
  peSRone_ += enq + en;
  return 1;
}

void PairHSCoulEwald::pairParticleParticleCheapEnergy_(const double &r2,
  const int &itype, const int &jtype, double * energy, double * force) {
  *energy = 0;
  const double sigij = sigij_[itype][jtype];
  if (r2 < sigij*sigij) {
    *energy = NUM_INF;
    peLJone_ += *energy;
  }
}

PairHSCoulEwald* PairHSCoulEwald::clone(Space* space) const {
  PairHSCoulEwald* p = new PairHSCoulEwald(*this);
  p->reconstruct(space);
  return p;
}

shared_ptr<PairHSCoulEwald> makePairHSCoulEwald(Space* space,
  const argtype &args) {
  return make_shared<PairHSCoulEwald>(space, args);
}

}  // namespace feasst
