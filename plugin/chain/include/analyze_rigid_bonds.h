
#ifndef FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_
#define FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

class AnalyzeRigidBonds : public AnalyzeUpdateOnly {
 public:
  AnalyzeRigidBonds(const argtype &args = argtype());
  void update(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    visitor_.compute(bond_, system.configuration());
    ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "bond check failure");
    visitor_.compute(angle_, system.configuration());
    ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "angle check failure");
  }

  std::string class_name() const override { return std::string("AnalyzeRigidBonds"); }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeRigidBonds>(istr); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
  }

  AnalyzeRigidBonds(std::istream& istr) : AnalyzeUpdateOnly(istr) {
    feasst_deserialize_version(istr); }

 private:
  BondVisitor visitor_;
  BondSquareWell bond_;
  AngleSquareWell angle_;
};

inline std::shared_ptr<AnalyzeRigidBonds> MakeAnalyzeRigidBonds(const argtype &args = argtype()) {
  return std::make_shared<AnalyzeRigidBonds>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_
