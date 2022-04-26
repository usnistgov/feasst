
#ifndef FEASST_STEPPERS_TUNE_PER_STATE_H_
#define FEASST_STEPPERS_TUNE_PER_STATE_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Each state and trial has a tunable value stored in this object.
  Every update, check if state has changed and replace with stored values.
  Also, Tune the stored values.
 */
class TunePerState : public Modify {
 public:
  /**
    args:
    - trials_per_tune: number of attempted trials per tune (default: 1e2).
    - stop_after_iteration: stop tuning when Criteria reaches this number
      of iterations. If -1, always tune (default: -1).
   */
  explicit TunePerState(argtype args = argtype());
  explicit TunePerState(argtype * args);

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  std::string write(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("TunePerState"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<TunePerState>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<TunePerState>(args); }
  explicit TunePerState(std::istream& istr);

 private:
  int trials_per_tune_;
  int stop_after_iteration_;
  std::vector<std::vector<double> > values_;
  std::vector<std::vector<int> > num_attempts_;
  std::vector<std::vector<int> > num_accepted_;
};

inline std::shared_ptr<TunePerState> MakeTunePerState(argtype args = argtype()) {
  return std::make_shared<TunePerState>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_TUNE_PER_STATE_H_
