
#ifndef FEASST_STEPPERS_CHECK_ENERGY_H_
#define FEASST_STEPPERS_CHECK_ENERGY_H_

#include "monte_carlo/include/modify.h"
#include "steppers/include/check.h"

namespace feasst {

/**
  Check if the running energy from criteria is equivalent, within tolerance,
  to an (unoptimized) calculation over the entire configuration.

  The origin of the energy deviation is different rounding at the level of
  numerical precision due to adding up the energy terms in different orders.
  For example, when an entire configuration energy is calculated, each pairwise
  energy term may be added in a different order than the sum of the contribution
  of each molecule.
  And each Monte Carlo trial may only calculate the energy change of one
  randomly chosen molecule.

  The energy check tolerance parameter is extensive, so it may be reasonable
  if the error happens when you increase the size.
  Thus, the tolerance parameter may need to be adjusted with system size.
  It is not very surprising to see a difference in the 6th or 7th decimal place
  in some cases.

  The size of the deviation can depend not only on the system size but the
  potential function and the number of trials between each check.
  For example, if the system undergoes a random walk deviation of size delta,
  it may deviate by the square root of the number of trials, where the size
  delta may be related to the number of interaction sites.
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
    Random * random,
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
