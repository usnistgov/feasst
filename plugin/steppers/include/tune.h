
#ifndef FEASST_STEPPERS_TUNE_H_
#define FEASST_STEPPERS_TUNE_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Each state and trial has a tunable value stored in this object.
  Every update, check if state has changed and replace with stored values.
  Also, Tune the stored values.
 */
class Tune : public Modify {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - trials_per_tune: number of attempted trials per tune (default: 1e3).
    - Stepper arguments.
   */
  explicit Tune(argtype args = argtype());
  explicit Tune(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  std::string write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("Tune"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<Tune>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<Tune>(args); }
  explicit Tune(std::istream& istr);

  //@}
 private:
  int trials_per_tune_;
  std::vector<double> values_;
  std::vector<int> num_attempts_;
  std::vector<int> num_accepted_;

  int min_num(const TrialFactory& trial_factory) const;
};

inline std::shared_ptr<Tune> MakeTune(argtype args = argtype()) {
  return std::make_shared<Tune>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_TUNE_H_
