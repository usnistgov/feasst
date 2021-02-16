
#ifndef FEASST_CHAIN_CHECK_RIGID_BONDS_H_
#define FEASST_CHAIN_CHECK_RIGID_BONDS_H_

#include "system/include/bond_visitor.h"
#include "system/include/bond_square_well.h"
#include "system/include/angle_square_well.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Assuming that all bonds are rigid, check that they are with "delta" (given as
  a bond parmeter) of the bond length.
 */
class CheckRigidBonds : public AnalyzeUpdateOnly {
 public:
  CheckRigidBonds(argtype args = argtype());
  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;
  std::string class_name() const override { return std::string("CheckRigidBonds"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<CheckRigidBonds>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit CheckRigidBonds(std::istream& istr);

 private:
  std::shared_ptr<BondVisitor> visitor_;
  BondSquareWell bond_;
  AngleSquareWell angle_;
};

inline std::shared_ptr<CheckRigidBonds> MakeCheckRigidBonds(
    argtype args = argtype()) {
  return std::make_shared<CheckRigidBonds>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_CHECK_RIGID_BONDS_H_
