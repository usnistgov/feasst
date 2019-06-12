
#ifndef FEASST_STEPPERS_TUNER_H_
#define FEASST_STEPPERS_TUNER_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Periodically adjust Trial parameters based on acceptance.
 */
class Tuner : public ModifyUpdateOnly {
 public:
  Tuner(const argtype &args = argtype()) : ModifyUpdateOnly(args) {}
  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    trial_factory->tune();
  }

  std::string class_name() const override { return std::string("Tuner"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<Tuner>(istr); }

  Tuner(std::istream& istr) : ModifyUpdateOnly(istr) {
    feasst_deserialize_version(istr); }

 private:
};

inline std::shared_ptr<Tuner> MakeTuner(const argtype &args = argtype()) {
  return std::make_shared<Tuner>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_TUNER_H_
