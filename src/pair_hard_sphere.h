/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
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
 * The diameter is determined by the sigma parameter described in Pair,
 * which is typically initialized by initData().
 */
class PairHardSphere : public Pair {
 public:
  /// Constructor
  PairHardSphere(Space* space);

  // Overloaded virtual function from pair.h
  // When pair parameters are initialized, automatically use sig2rCut.
  void initPairParam(const vector<double> eps,
    const vector<double> sig, const vector<double> sigref);

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
shared_ptr<PairHardSphere> makePairHardSphere(Space* space);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // PAIR_HARD_SPHERE_H_

