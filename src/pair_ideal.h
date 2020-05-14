/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_IDEAL_H_
#define PAIR_IDEAL_H_

#include <memory>
#include <vector>
#include "./pair.h"

namespace feasst {

/**
 * Ideal gas with no interactions.
 * PairIdeal still keeps track of neighbors for debugging purposes.
 */
class PairIdeal : public Pair {
 public:
  /// Constructor
  PairIdeal(Space* space, const argtype &args = argtype());

  // Construct from restart file
  PairIdeal(Space* space, const char* fileName);
  virtual ~PairIdeal() {}
  virtual PairIdeal* clone(Space* space) const;

 protected:
  // See comments of derived class from Pair
  void pairSiteSite_(const int &iSiteType, const int &jSiteType,
    double * energy, double * force, int * neighbor, const double &dx,
    const double &dy, const double &dz);

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairIdeal> makePairIdeal(Space* space,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // PAIR_IDEAL_H_

