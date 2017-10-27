/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_SQUAREWELL_H_
#define PAIR_SQUAREWELL_H_

#include <memory>
#include <vector>
#include "./pair.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

/**
 * Square well interactions.
 * The hard distance is determined by sigma,
 * the square well distance is determined by rCut,
 * and the depth of the well is epsilon.
 */
class PairSquareWell : public Pair {
 public:
  /// Constructor
  /// @param rCut interaciton cut off distance
  PairSquareWell(Space* space, const double rCut);

  /// Initialize hard sphere interactions
  /// between particle types itype and jtype.
  void initHardSphere(const int itype, const int jtype);

  // Overloaded virtual function from pair.h
  void multiPartEnerAtomCutInner(const double &r2, const int &itype,
                                 const int &jtype);

  // Construct from restart file
  PairSquareWell(Space* space, const char* fileName);
  virtual ~PairSquareWell() {}
  virtual PairSquareWell* clone(Space* space) const;

 protected:
  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairSquareWell> makePairSquareWell(Space* space, const double rCut);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_SQUAREWELL_H_

