
#ifndef FEASST_STEPPERS_CHECK_ENERGY_H_
#define FEASST_STEPPERS_CHECK_ENERGY_H_

#include "monte_carlo/include/modify.h"
#include "steppers/include/check.h"

namespace feasst {

/**
  Check that the running energy from criteria is equivalent, within tolerance,
  to an (unoptimized) calculation over the entire configuration.
 */
class CheckEnergy : public ModifyUpdateOnly {
 public:
  /**
    args:
    - tolerance: relative absolute difference between running energy
      and recomputed energy (default: 1e-10).
  */
  explicit CheckEnergy(argtype args = argtype());
  explicit CheckEnergy(argtype * args);

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  std::string class_name() const override { return std::string("CheckEnergy"); }

  // serialize
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CheckEnergy>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<CheckEnergy>(args); }
  explicit CheckEnergy(std::istream& istr);

 private:
  double tolerance_;
  std::shared_ptr<Analyze> check_;
};

inline std::shared_ptr<CheckEnergy> MakeCheckEnergy(
    argtype args = argtype()) {
  return std::make_shared<CheckEnergy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_ENERGY_H_
