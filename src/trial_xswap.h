/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_XSWAP_H_
#define TRIAL_XSWAP_H_

#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class TrialXSwap : public Trial {
 public:
  TrialXSwap();
  TrialXSwap(Space *space, Pair *pair, Criteria *criteria);
  TrialXSwap(const char* fileName, Space *space, Pair *pair,
             Criteria *criteria);
  ~TrialXSwap() {}
  TrialXSwap* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialXSwap* t = new TrialXSwap(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
  shared_ptr<TrialXSwap> cloneShrPtr(
    Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialXSwap, Trial>(
      cloneImpl(space, pair, criteria)));
  }
  void defaultConstruction();

  /// attempt random translation
  void attempt1();

 protected:
  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialXSwap> t = make_shared<TrialXSwap>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_XSWAP_H_

