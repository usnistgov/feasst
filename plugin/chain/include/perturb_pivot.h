
#ifndef FEASST_CHAIN_PERTURB_PIVOT_H_
#define FEASST_CHAIN_PERTURB_PIVOT_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

class PerturbPivot : public PerturbRotate {
 public:
  PerturbPivot(argtype args = argtype()) : PerturbPivot(&args) {
    FEASST_CHECK_ALL_USED(args);
  }
  PerturbPivot(argtype * args) : PerturbRotate(args) {
    class_name_ = "PerturbPivot";
  }

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot to the anchor.
  /// Dont rotate the particle positions.
  void move(const bool is_position_held, System * system, TrialSelect * select, Random * random) override;
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbPivot(std::istream& istr);
  virtual ~PerturbPivot() {}

 protected:
  void serialize_perturb_pivot_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbPivot> MakePerturbPivot(argtype args = argtype()) {
  return std::make_shared<PerturbPivot>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_PIVOT_H_
