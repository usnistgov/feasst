
#ifndef FEASST_CHAIN_PERTURB_ROTATE_COM_H_
#define FEASST_CHAIN_PERTURB_ROTATE_COM_H_

#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/**
  Uses the geometric center as the pivot point for the rotation.
 */
class PerturbRotateCOM : public PerturbRotate {
 public:
  explicit PerturbRotateCOM(argtype args = argtype());
  explicit PerturbRotateCOM(argtype * args);

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot to the geometric center of the mobile selection.
  void move(const bool is_position_held, System * system, TrialSelect * select,
    Random * random, Acceptance * acceptance) override;
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbRotateCOM(std::istream& istr);
  virtual ~PerturbRotateCOM() {}

 protected:
  void serialize_perturb_rotate_com_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbRotateCOM> MakePerturbRotateCOM(
    argtype args = argtype()) {
  return std::make_shared<PerturbRotateCOM>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_ROTATE_COM_H_
