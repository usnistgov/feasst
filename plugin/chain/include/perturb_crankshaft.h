
#ifndef FEASST_CHAIN_PERTURB_CRANKSHAFT_H_
#define FEASST_CHAIN_PERTURB_CRANKSHAFT_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

// HWH recenter particle position
class PerturbCrankshaft : public PerturbRotate {
 public:
  PerturbCrankshaft(const argtype& args = argtype()) : PerturbRotate(args) {
    class_name_ = "PerturbCrankshaft";
  }

  /// Set the pivot and axis of rotation by the ends of the selection.
  /// Select rotation angle randomly, bounded by tunable parameter.
  /// Dont rotate the particle positions.
  void move(System * system, TrialSelect * select, Random * random) override {
    const Position& pivot = select->mobile().site_positions()[0].front();
    axis_ = select->mobile().site_positions()[0].back();
    axis_.subtract(pivot);
    axis_.normalize();
    const double max_angle = tunable().value();
    const double angle = random->uniform_real(-max_angle, max_angle);
    rot_mat_.axis_angle(axis_, angle);
    PerturbRotate::move(pivot, rot_mat_, system, select,
      false // do not rotate particle positions
    );
  }

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
    const argtype& args = argtype()) {
  return std::make_shared<PerturbCrankshaft>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_CRANKSHAFT_H_
