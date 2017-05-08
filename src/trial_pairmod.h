/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef TRIAL_PAIRMOD_H_
#define TRIAL_PAIRMOD_H_

#include "./trial.h"

namespace feasst {

class TrialPairMod : public Trial {
 public:
  TrialPairMod();
  TrialPairMod(Space *space, Pair *pair, Criteria *criteria);
  TrialPairMod(const char* fileName, Space *space, Pair *pair,
               Criteria *criteria);
  ~TrialPairMod() {}
  TrialPairMod* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialPairMod* t = new TrialPairMod(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
  shared_ptr<TrialPairMod> cloneShrPtr(
    Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialPairMod, Trial>(
      cloneImpl(space, pair, criteria)));
  }

  /// attempt random translation
  void attempt1() {}

  /// ramp the chemical potential from value of min, located at
  // pair_->orderMin, to value of max for orderMax
  void initActivRamp(const double min, const double max) {
    lnzRampOn_ = 1; lnzMin_ = min; lnzMax_ = max;
  }

 protected:
  /// store the ramped chemical potential values
  int lnzRampOn_;
  double lnzMin_, lnzMax_;

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialPairMod> t = make_shared<TrialPairMod>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

}  // namespace feasst

#endif  // TRIAL_PAIRMOD_H_

