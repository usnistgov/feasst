/**
 * \file
 *
 * \brief trial add particles for Monte Carlo
 *
 */

#ifndef TRIAL_ADD_H_
#define TRIAL_ADD_H_

#include "./trial.h"

namespace feasst {

class Space;
class Pair;
class Criteria;

class TrialAdd : public Trial {
 public:
  explicit TrialAdd(const char* molType);
  TrialAdd(Space *space, Pair *pair, Criteria *criteria, const char* molType);
  TrialAdd(const char* fileName, Space *space, Pair *pair, Criteria *criteria);
  ~TrialAdd() {}
  TrialAdd* clone(Space* space, Pair *pair, Criteria *criteria) const;
  shared_ptr<TrialAdd> cloneShrPtr(Space* space, Pair* pair,
                                   Criteria* criteria) const;
  void defaultConstruction();
  void writeRestart(const char* fileName);

  /// attempt insertion
  void attempt1();

  /// return status of trial
  string printStat(const bool header = false);
  
  /// read only access to protected variables
  string molType() const { return molType_; }

 protected:
  string molType_;   //!< type of molecule to add

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(Space* space, Pair *pair,
                                      Criteria *criteria) const;
};

}  // namespace feasst

#endif  // TRIAL_ADD_H_

