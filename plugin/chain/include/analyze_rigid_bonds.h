
#ifndef FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_
#define FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_

#include "core/include/analyze.h"

namespace feasst {

class AnalyzeRigidBonds : public AnalyzeUpdateOnly {
 public:
  AnalyzeRigidBonds(const argtype &args = argtype()) : AnalyzeUpdateOnly(args) {}
  void update(const std::shared_ptr<Criteria> criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    visitor_.compute(bond_, system.configuration());
    ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "bond check failure");
    visitor_.compute(angle_, system.configuration());
    ASSERT(std::abs(visitor_.energy()) < NEAR_ZERO, "angle check failure");
  }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    auto model = std::make_shared<AnalyzeRigidBonds>();
    feasst_deserialize_version(istr);
    return model;
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

 private:
  const std::string class_name_ = "AnalyzeRigidBonds";
  BondVisitor visitor_;
  BondSquareWell bond_;
  AngleSquareWell angle_;
};

inline std::shared_ptr<AnalyzeRigidBonds> MakeAnalyzeRigidBonds(const argtype &args = argtype()) {
  return std::make_shared<AnalyzeRigidBonds>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_
