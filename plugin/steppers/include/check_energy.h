
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

  An alternative metric is the number of decimal places.
  In practice, this is implemented as follows for two energies U1 and U2.

  \f$\frac{|U_1-U_2|}{\mathrm{max}(|U_1|, |U_2|)} < 10^{-decimal\_places}\f$

  For this alternative metric, special consideration is made for small or near
  zero values of energy.
  In this special case, the check passes as long as

  \f$\mathrm{max}(|U_1|, |U_2|) < 10^{-decimal\_places}\f$

  Otherwise, the differences about zero to numerical precision would almost
  always fail the test.

  This class effectively functions as an Analyze because it does not change the
  System within the specified tolerance.
  However, the energy is recomputed and therefore the System is technically
  modified.
 */
class CheckEnergy : public ModifyUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - tolerance: absolute difference between running energy and recomputed
      energy (default: 1e-10).
    - decimal_places: If decimal places > 0 (default: 0), ignore the tolerance
      parameter and instead check if the difference between the two energies
      exceeds is within a the given number of decimal places.
    - Stepper arguments.
  */
  explicit CheckEnergy(argtype args = argtype());
  explicit CheckEnergy(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

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

  //@}
 private:
  double tolerance_;
  int decimal_places_;
  std::shared_ptr<Analyze> check_;

  bool is_within_tolerance_(const double u1, const double u2) const;
  std::string err_msg_() const;
};

inline std::shared_ptr<CheckEnergy> MakeCheckEnergy(
    argtype args = argtype()) {
  return std::make_shared<CheckEnergy>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_ENERGY_H_
