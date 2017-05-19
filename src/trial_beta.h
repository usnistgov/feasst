/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef TRIAL_BETA_H_
#define TRIAL_BETA_H_

#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class TrialBeta : public Trial {
 public:
  TrialBeta();
  TrialBeta(Space *space, Pair *pair, Criteria *criteria);
  TrialBeta(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialBeta() {}
  TrialBeta* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialBeta* t = new TrialBeta(*this); t->reconstruct(space, pair, criteria);
    return t; }
  shared_ptr<TrialBeta> cloneShrPtr(Space* space, Pair* pair,
                                    Criteria* criteria) const {
    return(std::static_pointer_cast<TrialBeta, Trial>
      (cloneImpl(space, pair, criteria))); }
  void defaultConstruction();

  /// attempt random translation
  void attempt1() {}

 protected:
  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const {
      shared_ptr<TrialBeta> t = make_shared<TrialBeta>(*this);
      t->reconstruct(space, pair, criteria);
      return t;}
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_BETA_H_

