
#ifndef FEASST_STEPPERS_CHECK_PHYSICALITY_H_
#define FEASST_STEPPERS_CHECK_PHYSICALITY_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically check that all sites are physical.
 */
class CheckPhysicality : public AnalyzeUpdateOnly {
 public:
  CheckPhysicality(const argtype &args = argtype()) : AnalyzeUpdateOnly(args) {}
  void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    for (int ic = 0; ic < system.num_configurations(); ++ic) {
      ASSERT(system.configuration(ic).are_all_sites_physical(),
        "all sites are not physical");
    }
  }

  std::string class_name() const override {
    return std::string("CheckPhysicality"); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CheckPhysicality>(istr); }

  CheckPhysicality(std::istream& istr) : AnalyzeUpdateOnly(istr) {
    feasst_deserialize_version(istr); }
};

inline std::shared_ptr<CheckPhysicality> MakeCheckPhysicality(
    const argtype &args = argtype()) {
  return std::make_shared<CheckPhysicality>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_CHECK_PHYSICALITY_H_
