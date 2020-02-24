
#ifndef FEASST_CHAIN_PERTURB_PIVOT_H_
#define FEASST_CHAIN_PERTURB_PIVOT_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

// HWH recenter particle position
class PerturbPivot : public PerturbRotate {
 public:
  PerturbPivot(const argtype& args = argtype()) : PerturbRotate(args) {
    class_name_ = "PerturbPivot";
  }

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot to the anchor.
  /// Dont rotate the particle positions.
  void move(System * system, TrialSelect * select, Random * random) override {
    const Position& pivot = select->anchor_position(0, 0, system);
    DEBUG("piv " << pivot.str());
    PerturbRotate::move(system, select, random, pivot, false);
    DEBUG(select->mobile().site_positions()[0][0].str());
  }
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbPivot(std::istream& istr);
  virtual ~PerturbPivot() {}

 protected:
  void serialize_perturb_pivot_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbPivot> MakePerturbPivot(const argtype& args = argtype()) {
  return std::make_shared<PerturbPivot>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_PIVOT_H_
