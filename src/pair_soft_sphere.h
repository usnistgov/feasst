/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_SOFT_SPHERE_H_
#define PAIR_SOFT_SPHERE_H_

#include <memory>
#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Soft sphere pair-wise interaction, \f$ U=eps(sig/r)^n \f$
 */
class PairSoftSphere : public Pair {
 public:
  /// Constructor
  /// @param rCut interaciton cut off distance
  PairSoftSphere(Space* space, const double rCut);

  /// Initialize the potential parameter exponent (default below).
  void initExponent(const int n = 12) { n_ = n; }

  /// Return the exponential parameter.
  double exponent() const { return n_; }

  /// Return the reduced second virial coefficient,
  /// b2~ = b2(beta eps)^(-3/n==12) as a reference
  double b2reduced();

  // Overloaded virtual function from pair.h
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  // Write restart file.
  virtual void writeRestart(const char* fileName);

  // Construct from restart file
  PairSoftSphere(Space* space, const char* fileName);
  virtual ~PairSoftSphere() {}
  virtual PairSoftSphere* clone(Space* space) const;

 protected:
  double n_;    // potential parameter for power order

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairSoftSphere> makePairSoftSphere(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_SOFT_SPHERE_H_

