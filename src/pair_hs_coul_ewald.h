/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_HS_COUL_EWALD_H_
#define PAIR_HS_COUL_EWALD_H_

#include <memory>
#include <vector>
#include "./pair_lj_coul_ewald.h"

namespace feasst {

/**
 * Hard sphere pair-wise interaction.
 * The diameter is determined by the sigma parameter described in Pair,
 * which is typically initialized by initData().
 */
class PairHSCoulEwald : public PairLJCoulEwald {
 public:
  /// Constructor
  PairHSCoulEwald(Space* space, const argtype &args = argtype());

  // Construct from restart file
  PairHSCoulEwald(Space* space, const char* fileName);
  virtual ~PairHSCoulEwald() {}
  virtual PairHSCoulEwald* clone(Space* space) const;

  // Overloaded virtual function from pair.h
  int multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                const int &jtype);

  // Overloaded virtual in order to use LJ of Oxygen when cheap energy enabled
  void pairParticleParticleCheapEnergy_(const double &r2, const int &itype,
    const int &jtype, double * energy, double * force);

 protected:
  // defaults in constructor
  void defaultConstruction_();

  // See comments of derived class from Pair
  void pairSiteSite_(const int &iSiteType, const int &jSiteType,
    double * energy, double * force, int * neighbor, const double &dx,
    const double &dy, const double &dz);
};

/// Factory method
shared_ptr<PairHSCoulEwald> makePairHSCoulEwald(Space* space,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // PAIR_HS_COUL_EWALD_H_
