/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST’s creation of the data/software.
 */

#ifndef PAIR_IDEAL_H_
#define PAIR_IDEAL_H_

#include <memory>
#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Ideal gas with no interactions.
 * PairIdeal still keeps track of neighbors for debugging purposes.
 */
class PairIdeal : public Pair {
 public:
  /// Constructor
  /// @param rCut HWH:depreciated. This variable has no use in this class.
  PairIdeal(Space* space, const double rCut);

  // Overloaded virtual function from pair.h
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  // Construct from restart file
  PairIdeal(Space* space, const char* fileName);
  virtual ~PairIdeal() {}
  virtual PairIdeal* clone(Space* space) const;

 protected:
  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairIdeal> makePairIdeal(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_IDEAL_H_

