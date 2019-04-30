
#ifndef FEASST_CORE_TUNER_H_
#define FEASST_CORE_TUNER_H_

#include "core/include/modify.h"

namespace feasst {

/**
 */
class Tuner : public ModifyUpdateOnly {
 public:
  Tuner(const argtype &args = argtype()) : ModifyUpdateOnly(args) {}
  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    trial_factory->tune();
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    feasst_deserialize_version(istr);
    auto modify = std::make_shared<Tuner>();
    return modify;
  }

 private:
  const std::string class_name_ = "Tuner";
};

inline std::shared_ptr<Tuner> MakeTuner(const argtype &args = argtype()) {
  return std::make_shared<Tuner>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_TUNER_H_
