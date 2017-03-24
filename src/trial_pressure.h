/**
 * \file
 *
 * \brief trial moves for Monte Carlo
 *
 */

#ifndef TRIAL_PRESSURE_H_
#define TRIAL_PRESSURE_H_

#include "./trial.h"

class TrialPressure : public Trial {
 public:
  explicit TrialPressure(const char* variable);
  TrialPressure(Space *space, Pair *pair, Criteria *criteria,
                const char* variable);
  TrialPressure(const char* fileName, Space *space, Pair *pair,
                Criteria *criteria);
  ~TrialPressure() {}
  TrialPressure* clone(Space* space, Pair* pair, Criteria* criteria) const {
    TrialPressure* t = new TrialPressure(*this);
    t->reconstruct(space, pair, criteria); return t;
  }
  shared_ptr<TrialPressure> cloneShrPtr(
    Space* space, Pair* pair, Criteria* criteria) const {
    return(std::static_pointer_cast<TrialPressure, Trial>(
      cloneImpl(space, pair, criteria)));
  }

  /// attempt random translation
  void attempt1() {}

 protected:
  /// pressure variable (e.g., pressure, or lnpres [ln(pressure)])
  string variable_;

  // clone design pattern
  virtual shared_ptr<Trial> cloneImpl(
    Space* space, Pair *pair, Criteria *criteria) const {
    shared_ptr<TrialPressure> t = make_shared<TrialPressure>(*this);
    t->reconstruct(space, pair, criteria);
    return t;
  }
};

#endif  // TRIAL_PRESSURE_H_

