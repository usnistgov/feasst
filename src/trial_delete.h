/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_DELETE_H_
#define TRIAL_DELETE_H_

#include <memory>
#include "./trial.h"

class Space;
class Pair;
class Criteria;

class TrialDelete : public Trial {
 public:
  TrialDelete();
  explicit TrialDelete(const char* molType);
  TrialDelete(Space *space, Pair *pair, Criteria *criteria);
  TrialDelete(Space *space, Pair *pair, Criteria *criteria, const char* molType);
  TrialDelete(const char* fileName, Space *space, Pair *pair,
              Criteria *criteria);
  ~TrialDelete() {}
  TrialDelete* clone(Space* space, Pair *pair, Criteria *criteria) const {
    TrialDelete* t = new TrialDelete(*this);
    t->reconstruct(space, pair, criteria); return t; }
  shared_ptr<TrialDelete> cloneShrPtr
    (Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialDelete, Trial>
    (cloneImpl(space, pair, criteria))); }

  // default constructor
  void defaultConstruction();
  void writeRestart(const char* fileName);

  /// attempt deletion
  void attempt1();

 protected:
  string molType_;   //!< type of molecule to delete

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl
    (Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialDelete> t = make_shared<TrialDelete>(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
};

#endif  // TRIAL_DELETE_H_

