
#ifndef FEASST_STEPPERS_CHECK_PHYSICALITY_H_
#define FEASST_STEPPERS_CHECK_PHYSICALITY_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically check that all sites are physical.
 */
class CheckPhysicality : public AnalyzeUpdateOnly {
 public:
  CheckPhysicality(argtype args = argtype());
  CheckPhysicality(argtype * args);
  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;
  std::string class_name() const override {
    return std::string("CheckPhysicality"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CheckPhysicality>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<CheckPhysicality>(args); }
  explicit CheckPhysicality(std::istream& istr);
};

inline std::shared_ptr<CheckPhysicality> MakeCheckPhysicality(
    argtype args = argtype()) {
  return std::make_shared<CheckPhysicality>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_PHYSICALITY_H_
