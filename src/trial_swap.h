/**
 * \file
 *
 * \brief trial add particles for Monte Carlo
 *
 */

#ifndef TRIAL_SWAP_H_
#define TRIAL_SWAP_H_

#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class Space;
class Pair;
class Criteria;

class TrialSwap : public Trial {
 public:
  explicit TrialSwap(const char* molTypeA, const char* molTypeB);
  TrialSwap(Space *space, Pair *pair, Criteria *criteria, const char* molTypeA, const char* molTypeB);
  TrialSwap(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialSwap() {}
  TrialSwap* clone(Space* space, Pair *pair, Criteria *criteria) const;
  shared_ptr<TrialSwap> cloneShrPtr(Space* space, Pair* pair,
                                   Criteria* criteria) const;
  void defaultConstruction();
  void writeRestart(const char* fileName);

  void attempt1();

  /// return status of trial
  string printStat(const bool header = false);

  /// read only access to protected variables
  string molTypeA() const { return molTypeA_; }
  string molTypeB() const { return molTypeB_; }
  Accumulator nA() const { return nA_; };
  Accumulator nB() const { return nB_; };

 protected:
  string molTypeA_;
  string molTypeB_;
  Accumulator nA_, nB_;

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

