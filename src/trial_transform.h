/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_TRANSFORM_H_
#define TRIAL_TRANSFORM_H_

#include "./trial.h"

class TrialTransform : public Trial {
 public:
  explicit TrialTransform(const char* transType);
  TrialTransform(Space *space, Pair *pair, Criteria *criteria,
                 const char* transType);
  TrialTransform(const char* fileName, Space *space, Pair *pair,
                 Criteria *criteria);
  ~TrialTransform() {}
  TrialTransform* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialTransform* t = new TrialTransform(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
  shared_ptr<TrialTransform> cloneShrPtr(
    Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialTransform, Trial>(
      cloneImpl(space, pair, criteria)));
  }
  void defaultConstruction();
  void writeRestart(const char* fileName);

  /// attempt random translation
  void attempt1();

  // tune parameters (e.g., based on acceptance)
  void tuneParameters();
  double targAcceptPer;      //!< target acceptance percentage

  /// return status of trial
  string printStat(const bool header = false);
  
  // functions for read-only access of private data-members
  string transType() const { return transType_; }

 protected:
  string transType_;  //!< type of transformation

  void scaleAttempt_(const double factor);
  
  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialTransform> t = make_shared<TrialTransform>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

#endif  // TRIAL_TRANSFORM_H_

