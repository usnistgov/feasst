/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef TRIAL_CLUSTER_H_
#define TRIAL_CLUSTER_H_

#include "./trial.h"

namespace feasst {

class TrialCluster : public Trial {
 public:
  explicit TrialCluster(const char* transType);
  TrialCluster(Space *space, Pair *pair, Criteria *criteria,
               const char* transType);
  TrialCluster(const char* fileName, Space *space, Pair *pair,
               Criteria *criteria);
  ~TrialCluster() {}
  TrialCluster* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialCluster* t = new TrialCluster(*this);
    t->reconstruct(space, pair, criteria); return t; }
  shared_ptr<TrialCluster> cloneShrPtr
    (Space* space, Pair* pair, Criteria* criteria) const {
      return(std::static_pointer_cast<TrialCluster, Trial>
      (cloneImpl(space, pair, criteria))); }

  /// attempt random translation
  void attempt1() {}

  // tune parameters (e.g., based on acceptance)
  void tuneParameters() {}
  double targAcceptPer;      //!< target acceptance percentage
  double clusterCut;      //!< cut-off for cluster definition

  // functions for read-only access of private data-members
  string transType() const { return transType_; }

 protected:
  string transType_;  //!< type of transformation

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const {
      shared_ptr<TrialCluster> t = make_shared<TrialCluster>(*this);
      t->reconstruct(space, pair, criteria);
      return t;
  }
};

}  // namespace feasst

#endif  // TRIAL_CLUSTER_H_

