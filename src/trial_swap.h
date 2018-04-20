/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef TRIAL_SWAP_H_
#define TRIAL_SWAP_H_

#include <memory>
#include <string>
#include "./trial.h"

namespace feasst {

class Space;
class Pair;
class Criteria;

/**
 * Attempt to randomly swap the type of a molecule.
 * As currently implemented, the position of the first particle in the molecule
 * during swap is where the first particle in the other molecule is placed.
 * Thus, for multi-site particles, the orientation is random.
 */
class TrialSwap : public Trial {
 public:
  /// Constructor with particle types which are allowed to swap.
  TrialSwap(Pair *pair, Criteria *criteria, const char* molTypeA,
            const char* molTypeB);

  /// This constructor is not often used, but its purpose is to initialize trial
  /// for interface before using reconstruct to set object pointers.
  TrialSwap(const char* molTypeA, const char* molTypeB);

  /// Return the first type of particle allowed to swap.
  string molTypeA() const { return molTypeA_; }

  /// Return the Accumulator for the number of "A" particles.
  Accumulator nA() const { return nA_; }

  /// Return the second type of particle allowed to swap.
  string molTypeB() const { return molTypeB_; }

  /// Return the Accumulator for the number of "B" particles.
  Accumulator nB() const { return nB_; }

  /// return status of trial
  string printStat(const bool header = false);

  /// Write restart file.
  void writeRestart(const char* fileName);

  /// Construct from restart file.
  TrialSwap(const char* fileName, Pair *pair, Criteria *criteria);
  ~TrialSwap() {}
  TrialSwap* clone(Pair* pair, Criteria* criteria) const {
    TrialSwap* t = new TrialSwap(*this);
    t->reconstruct(pair, criteria); return t;
  }
  shared_ptr<TrialSwap> cloneShrPtr(
    Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialSwap, Trial>(
      cloneImpl(pair, criteria)));
  }

 protected:
  string molTypeA_;
  string molTypeB_;
  Accumulator nA_, nB_;

  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialSwap> t = make_shared<TrialSwap>(*this);
    t->reconstruct(pair, criteria);
    return t;
  }
};

/// Factory method
shared_ptr<TrialSwap> makeTrialSwap(Pair *pair,
  Criteria *criteria, const char* molTypeA, const char* molTypeB);

/// Factory method
shared_ptr<TrialSwap> makeTrialSwap(const char* molTypeA, const char* molTypeB);

class MC;
void swapTrial(MC *mc, const char* molTypeA, const char* molTypeB);
void swapTrial(shared_ptr<MC> mc, const char* molTypeA, const char* molTypeB);

}  // namespace feasst

#endif  // TRIAL_SWAP_H_

