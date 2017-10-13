/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_HARD_SPHERE_H_
#define PAIR_HARD_SPHERE_H_

#include <memory>
#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Hard sphere pair-wise interaction.
 * The diameter is determined by the sigma parameter described in Pair.
 * Use the convenience function pair.sig2rCut() to set the rCut's to sigma.
 */
class PairHardSphere : public Pair {
 public:
  /// Constructor
  /// @param rCut interaciton cut off distance should be equal to sigma
  PairHardSphere(Space* space, const double rCut);

  // Overloaded virtual function from pair.h
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  // Construct from restart file
  PairHardSphere(Space* space, const char* fileName);
  virtual ~PairHardSphere() {}
  virtual PairHardSphere* clone(Space* space) const;

 protected:
  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairHardSphere> makePairHardSphere(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HARD_SPHERE_H_

