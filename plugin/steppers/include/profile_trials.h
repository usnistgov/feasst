
#ifndef FEASST_STEPPERS_PROFILE_TRIALS_H_
#define FEASST_STEPPERS_PROFILE_TRIALS_H_

#include <ctime>
#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Output the percentage of CPU time spent on each Trial.
  Consider using trials_per_update != 1,
  because profiling may be time consuming. 
 */
class ProfileTrials : public Analyze {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit ProfileTrials(argtype args = argtype());
  explicit ProfileTrials(argtype * args);

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

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  const std::vector<Accumulator>& profile() const { return profile_; }

  // serialize
  std::string class_name() const override { return std::string("ProfileTrials"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<ProfileTrials>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<ProfileTrials>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ProfileTrials(std::istream& istr);

  //@}
 private:
  std::vector<Accumulator> profile_;

  // temporary and not to be serialized
  clock_t previous_clock_;
  bool is_previous_ = false;
};

inline std::shared_ptr<ProfileTrials> MakeProfileTrials(
    argtype args = argtype()) {
  return std::make_shared<ProfileTrials>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_PROFILE_TRIALS_H_
