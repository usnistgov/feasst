/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef TRIAL_GCA_H_
#define TRIAL_GCA_H_

#include "./trial.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

class TrialGCA : public Trial {
 public:
  TrialGCA();
  TrialGCA(Space *space, Pair *pair, Criteria *criteria);
  TrialGCA(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialGCA() {}
  TrialGCA* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialGCA* t = new TrialGCA(*this); t->reconstruct(space, pair, criteria);
    return t; }
  shared_ptr<TrialGCA> cloneShrPtr(
    Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialGCA, Trial>(
      cloneImpl(space, pair, criteria)));
  }

  /// attempt random translation
  void attempt1() {}

  // tune parameters (e.g., based on acceptance)
  void tuneParameters() {}
  double targAcceptPer;      //!< target acceptance percentage

 protected:
  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialGCA> t = make_shared<TrialGCA>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // TRIAL_GCA_H_

