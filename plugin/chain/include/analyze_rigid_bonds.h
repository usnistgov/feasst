
#ifndef FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_
#define FEASST_CHAIN_ANALYZE_RIGID_BONDS_H_

#include "system/include/bond_visitor.h"
#include "system/include/bond_square_well.h"
#include "system/include/angle_square_well.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

class AnalyzeRigidBonds : public AnalyzeUpdateOnly {
 public:
  AnalyzeRigidBonds(const argtype &args = argtype());
  void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;
  std::string class_name() const override { return std::string("AnalyzeRigidBonds"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeRigidBonds>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit AnalyzeRigidBonds(std::istream& istr);

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
