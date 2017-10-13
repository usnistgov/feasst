/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_HARD_CIRCLE_H_
#define PAIR_HARD_CIRCLE_H_

#include <memory>
#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Hard circles with implicit depletant model (Asakura-Oosawa) based on
 * overlap of excluded area.
 * \f$ U = U_{Hard} + \frac{\Delta A_{ex}}{\pi R_g^2} \phi k_B T\f$
 * The hard particle diameter is set to 1. by default (dCircle_).
 */
class PairHardCircle : public Pair {
 public:
  /// Constructor
  /// @param rCut interaciton cut off distance
  PairHardCircle(Space* space, const double rCut);

  /// Initialize the radius of gyration of the depletant, \f$R_g\f$.
  void initRDep(const double rDep = 0.) { rDep_ = rDep; }

  // Overloaded virtual function from pair.h
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  // Write restart file.
  virtual void writeRestart(const char* fileName);

  // Construct from restart file
  PairHardCircle(Space* space, const char* fileName);
  virtual ~PairHardCircle() {}
  virtual PairHardCircle* clone(Space* space) const;

 protected:
  double dCircle_;               //!< diameter of hard circle
  double rDep_;               //!< radius of depletant

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairHardCircle> makePairHardCircle(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HARD_CIRCLE_H_

