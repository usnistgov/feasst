
#ifndef FEASST_STEPPERS_TUNER_H_
#define FEASST_STEPPERS_TUNER_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Periodically adjust Trial parameters based on acceptance.
 */
class Tuner : public ModifyUpdateOnly {
 public:
  explicit Tuner(argtype args = argtype());

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override {
    trial_factory->tune(); }

  // serialize
  std::string class_name() const override { return std::string("Tuner"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<Tuner>(istr); }
  explicit Tuner(std::istream& istr);
};

inline std::shared_ptr<Tuner> MakeTuner(argtype args = argtype()) {
  return std::make_shared<Tuner>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_TUNER_H_
