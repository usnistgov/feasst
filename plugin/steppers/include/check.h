
#ifndef FEASST_STEPPERS_CHECK_H_
#define FEASST_STEPPERS_CHECK_H_

#include <memory>
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Run periodic checks on each class.
 */
class Check : public AnalyzeUpdateOnly {
 public:
  Check(const argtype &args = argtype()) : AnalyzeUpdateOnly(args) {}
  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;
  std::string class_name() const override { return std::string("Check"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Check>(istr); }
  explicit Check(std::istream& istr);
};

inline std::shared_ptr<Check> MakeCheck(const argtype &args = argtype()) {
  return std::make_shared<Check>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_H_
