/**
 * This file is a stub or placeholder for an experimental class that is not part of this release.
 */

#ifndef TRIAL_CONFIGBIAS_H_
#define TRIAL_CONFIGBIAS_H_

#include "./trial.h"

namespace feasst {

class Space;
class Pair;
class Criteria;
class MC;

class TrialConfigBias : public Trial {
 public:
  explicit TrialConfigBias(const char* molType);
  TrialConfigBias(Space *space, Pair *pair, Criteria *criteria,
                  const char* molType);
  TrialConfigBias(const char* fileName, Space *space, Pair *pair,
                  Criteria *criteria);
  ~TrialConfigBias() {}
  TrialConfigBias* clone(Space* space, Pair *pair, Criteria *criteria) const;
  shared_ptr<TrialConfigBias> cloneShrPtr
    (Space* space, Pair* pair, Criteria* criteria) const;

  /// attempt insertion
  void attempt1() {}

  /// flag on how to handle insert/delete trials
  int insDelFlag;
  
  /// initialize trials in given MC class based on reading input file
  void file2MC(const char* fileName, MC* mc) {}
 
 protected:
  string molType_;              //!< type of molecule to add
  
  /// clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const;
};

}  // namespace feasst

#endif  // TRIAL_CONFIGBIAS_H_

