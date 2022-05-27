
#ifndef FEASST_CHAIN_PERTURB_CRANKSHAFT_SMALL_H_
#define FEASST_CHAIN_PERTURB_CRANKSHAFT_SMALL_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/// A more efficient version of PerturbCrankshaft for small molecules.
class PerturbCrankshaftSmall : public PerturbRotate {
 public:
  PerturbCrankshaftSmall(argtype args = argtype()) : PerturbCrankshaftSmall(&args) {
    FEASST_CHECK_ALL_USED(args);
  }
  PerturbCrankshaftSmall(argtype * args) : PerturbRotate(args) {
    class_name_ = "PerturbCrankshaftSmall";
  }

  /// Set the pivot and axis of rotation by the anchors.
  /// Select rotation angle randomly, bounded by tunable parameter.
  void move(const bool is_position_held, System * system, TrialSelect * select, Random * random) override;
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbCrankshaftSmall(std::istream& istr);
  virtual ~PerturbCrankshaftSmall() {}

 protected:
  void serialize_perturb_crankshaft_(std::ostream& ostr) const;

 private:
  // temporary
  Position axis_;
  RotationMatrix rot_mat_;
};

inline std::shared_ptr<PerturbCrankshaftSmall> MakePerturbCrankshaftSmall(
    argtype args = argtype()) {
  return std::make_shared<PerturbCrankshaftSmall>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_CRANKSHAFT_SMALL_H_
