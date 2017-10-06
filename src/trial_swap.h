#ifndef TRIAL_SWAP_H_
#define TRIAL_SWAP_H_

#include <memory>
#include <string>
#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

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
  TrialSwap(Space *space, Pair *pair, Criteria *criteria, const char* molTypeA,
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
  TrialSwap(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialSwap() {}
  TrialSwap* clone(Space* space, Pair *pair, Criteria *criteria) const;
  shared_ptr<TrialSwap> cloneShrPtr(Space* space, Pair* pair,
                                   Criteria* criteria) const;

 protected:
  string molTypeA_;
  string molTypeB_;
  Accumulator nA_, nB_;

  void attempt1_();

  void defaultConstruction_();

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(Space* space, Pair *pair,
                                      Criteria *criteria) const;
};

class MC;
void swapTrial(MC *mc, const char* molTypeA, const char* molTypeB);
void swapTrial(shared_ptr<MC> mc, const char* molTypeA, const char* molTypeB);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_SWAP_H_

