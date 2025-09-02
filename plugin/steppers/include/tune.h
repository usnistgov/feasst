
#ifndef FEASST_STEPPERS_TUNE_H_
#define FEASST_STEPPERS_TUNE_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Each trial has a tunable value stored in this object.
  Every update, check if state has changed and replace with stored values.
  Also, Tune the stored values.
 */
class Tune : public Modify {
 public:
  //@{
  /** @name Arguments
    - trials_per_tune: number of attempted trials per tune (default: 1e3).
    - Stepper arguments.
   */
  explicit Tune(argtype args = argtype());
  explicit Tune(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(MonteCarlo * mc) override;
  std::string write(MonteCarlo * mc) override;

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

  std::vector<double> * get_values_() { return &((*data_.get_dble_2D())[0]); }
  const std::vector<double>& values_() const { return data_.dble_2D()[0]; }
  std::vector<int> * get_num_attempts_() { return &((*data_.get_int_2D())[0]); }
  const std::vector<int>& num_attempts_() const { return data_.int_2D()[0]; }
  std::vector<int> * get_num_accepted_() { return &((*data_.get_int_2D())[1]); }
  const std::vector<int>& num_accepted_() const { return data_.int_2D()[1]; }
  int min_num(const TrialFactory& trial_factory) const;
};

inline std::shared_ptr<Tune> MakeTune(argtype args = argtype()) {
  return std::make_shared<Tune>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_TUNE_H_
