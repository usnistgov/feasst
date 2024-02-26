
#ifndef FEASST_CHAIN_PERTURB_CRANKSHAFT_H_
#define FEASST_CHAIN_PERTURB_CRANKSHAFT_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

class PerturbCrankshaft : public PerturbRotate {
 public:
  PerturbCrankshaft(argtype args = argtype()) : PerturbCrankshaft(&args) {
    FEASST_CHECK_ALL_USED(args);
  }
  PerturbCrankshaft(argtype * args) : PerturbRotate(args) {
    class_name_ = "PerturbCrankshaft";
  }

  /// Set the pivot and axis of rotation by the ends of the selection.
  /// Select rotation angle randomly, bounded by tunable parameter.
  /// Dont rotate the particle positions.
  void move(const bool is_position_held, System * system, TrialSelect * select,
    Random * random, Acceptance * acceptance) override;
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbCrankshaft(std::istream& istr);
  virtual ~PerturbCrankshaft() {}

 protected:
  void serialize_perturb_crankshaft_(std::ostream& ostr) const;

 private:
  // temporary
  Position axis_;
  RotationMatrix rot_mat_;
};

inline std::shared_ptr<PerturbCrankshaft> MakePerturbCrankshaft(
    argtype args = argtype()) {
  return std::make_shared<PerturbCrankshaft>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_CRANKSHAFT_H_
