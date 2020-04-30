
#ifndef FEASST_CLUSTER_PERTURB_MOVE_AVB_H_
#define FEASST_CLUSTER_PERTURB_MOVE_AVB_H_

#include "monte_carlo/include/perturb_move.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

/**
  Move Select to inside or outside the AVB of anchor.
 */
class PerturbMoveAVB : public PerturbMove {
 public:
  /**
    args:
    - inside: true if selecting in the AV, otherwise out (default: true).
      Not implemented for grand_canonical.
   */
  PerturbMoveAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());

  void move(System * system,
            TrialSelect * select,
            Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbMoveAVB(std::istream& istr);
  virtual ~PerturbMoveAVB() {}

 private:
  std::shared_ptr<NeighborCriteria> neighbor_criteria_;
  PerturbRotate rotate_;
  bool inside_;

  // temporary
  PerturbTranslate translate_;
  Position displace_;
};

inline std::shared_ptr<PerturbMoveAVB> MakePerturbMoveAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype()) {
  return std::make_shared<PerturbMoveAVB>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_PERTURB_MOVE_AVB_H_
