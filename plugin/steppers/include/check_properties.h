
#ifndef FEASST_STEPPERS_CHECK_PROPERTIES_H_
#define FEASST_STEPPERS_CHECK_PROPERTIES_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Check equivalence of stored properties and updated properties.
 */
class CheckProperties : public ModifyUpdateOnly {
 public:
  /**
    args:
    - tolerance: acceptable deviation (default: 1e-15).
   */
  explicit CheckProperties(argtype args = argtype());
  explicit CheckProperties(argtype * args);

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override {
    return std::string("CheckProperties"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CheckProperties>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<CheckProperties>(args); }
  explicit CheckProperties(std::istream& istr);

 private:
  double tolerance_;
};

inline std::shared_ptr<CheckProperties> MakeCheckProperties(
    argtype args = argtype()) {
  return std::make_shared<CheckProperties>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_PROPERTIES_H_
