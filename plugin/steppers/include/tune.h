
#ifndef FEASST_STEPPERS_TUNE_H_
#define FEASST_STEPPERS_TUNE_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Periodically adjust Trial parameters based on acceptance.
 */
class Tune : public ModifyUpdateOnly {
 public:
  explicit Tune(argtype args = argtype());
  explicit Tune(argtype * args);

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override {
    trial_factory->tune(); }

  // serialize
  std::string class_name() const override { return std::string("Tune"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<Tune>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<Tune>(args); }
  explicit Tune(std::istream& istr);
};

inline std::shared_ptr<Tune> MakeTune(argtype args = argtype()) {
  return std::make_shared<Tune>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_TUNE_H_
