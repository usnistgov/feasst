
#ifndef FEASST_CHAIN_PERTURB_ROTATE_COM_H_
#define FEASST_CHAIN_PERTURB_ROTATE_COM_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

// HWH recenter particle position
class PerturbRotateCOM : public PerturbRotate {
 public:
  PerturbRotateCOM(const argtype& args = argtype()) : PerturbRotate(args) {
    class_name_ = "PerturbRotateCOM";
  }

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot to the geometric center of the mobile selection.
  /// Rotate the particle positions.
  void move(System * system, TrialSelect * select, Random * random) override {
    const Position pivot = select->mobile().geometric_center();
    DEBUG("piv " << pivot.str());
    PerturbRotate::move(system, select, random, pivot, true);
  }
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbRotateCOM(std::istream& istr);
  virtual ~PerturbRotateCOM() {}

 protected:
  void serialize_perturb_rotate_com_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbRotateCOM> MakePerturbRotateCOM(
    const argtype& args = argtype()) {
  return std::make_shared<PerturbRotateCOM>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_ROTATE_COM_H_
