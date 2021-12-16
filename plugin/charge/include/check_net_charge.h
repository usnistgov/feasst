
#ifndef FEASST_CHARGE_CHECK_NET_CHARGE_H_
#define FEASST_CHARGE_CHECK_NET_CHARGE_H_

#include "monte_carlo/include/analyze.h"
#include "charge/include/ewald.h"

namespace feasst {

/**
  Periodically check if charge of the system is within an acceptable range.
 */
class CheckNetCharge : public AnalyzeUpdateOnly {
 public:
  /**
    args:
    - minimum: minimum acceptable charge (default: 0).
    - maximum: maximum acceptable charge (default: 0).
   */
  CheckNetCharge(argtype args = argtype());
  void update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) override;
  std::string class_name() const override {
    return std::string("CheckNetCharge"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CheckNetCharge>(istr); }
  explicit CheckNetCharge(std::istream& istr);
  virtual ~CheckNetCharge() {}

 protected:
  Ewald ewald_;
  double minimum_;
  double maximum_;
};

inline std::shared_ptr<CheckNetCharge> MakeCheckNetCharge(
    argtype args = argtype()) {
  return std::make_shared<CheckNetCharge>(args);
}

}  // namespace feasst

#endif  // FEASST_CHARGE_CHECK_NET_CHARGE_H_
