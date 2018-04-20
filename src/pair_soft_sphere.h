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

namespace feasst {

/**
 * Soft sphere pair-wise interaction, \f$ U=eps(sig/r)^n \f$
 */
class PairSoftSphere : public Pair {
 public:
  /// Constructor
  PairSoftSphere(Space* space, const argtype &args = argtype());

  /// Initialize the potential parameter exponent (default below).
  void initExponent(const int n = 12) { n_ = n; }

  /// Return the exponential parameter.
  double exponent() const { return n_; }

  /// Return the reduced second virial coefficient,
  /// b2~ = b2(beta eps)^(-3/n==12) as a reference
  double b2reduced();

  // Write restart file.
  virtual void writeRestart(const char* fileName);

  // Construct from restart file
  PairSoftSphere(Space* space, const char* fileName);
  virtual ~PairSoftSphere() {}
  virtual PairSoftSphere* clone(Space* space) const;

 protected:
  double n_;    // potential parameter for power order

  // See comments of derived class from Pair
  void pairSiteSite_(const int &iSiteType, const int &jSiteType,
    double * energy, double * force, int * neighbor, const double &dx,
    const double &dy, const double &dz);


  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairSoftSphere> makePairSoftSphere(Space* space,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // PAIR_SOFT_SPHERE_H_

