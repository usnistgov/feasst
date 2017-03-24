/**
 * \file
 *
 * \brief
 *
 */

#ifndef TRIAL_GROW_H_
#define TRIAL_GROW_H_

#include "./trial.h"
#include <memory>
#include <string>
#include <vector>

class TrialGrow : public Trial {
 public:
  TrialGrow(const char* molType, const int nStages);
  TrialGrow(Space *space, Pair *pair, Criteria *criteria, const char* molType,
            const int nStages);
  TrialGrow(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialGrow() {}
  TrialGrow* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialGrow* t = new TrialGrow(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
  shared_ptr<TrialGrow> cloneShrPtr(Space* space, Pair* pair, Criteria* criteria
    ) const { return(std::static_pointer_cast<TrialGrow, Trial>(
      cloneImpl(space, pair, criteria)));
  }

  /// randomly attempt grow or shrink move with equal probability
  void attempt1() {}

  // functions for read-only access of private data-members
  string molType() const { return molType_; }
  int nStages() const { return static_cast<int>(stageScale_.size()); }

 protected:
  string molType_;              //!< type of molecule to grow
  int currentStage_;            //!< current stage of growth
  vector<double> stageScale_;   //!< scale of stages

  /// for each stage, distances of atoms from first atom for original molecule
  //  (100% scale)
  vector<vector<double> > bondLengths_;

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(Space* space, Pair *pair,
    Criteria *criteria) const {
    shared_ptr<TrialGrow> t = make_shared<TrialGrow>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

#endif  // TRIAL_GROW_H_

