
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
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
   */
  explicit PerturbMoveAVB(argtype args = argtype());
  explicit PerturbMoveAVB(argtype * args);

  void move(const bool is_position_held,
            System * system,
            TrialSelect * select,
            Random * random,
            Acceptance * acceptance) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbMoveAVB(std::istream& istr);
  virtual ~PerturbMoveAVB() {}

 private:
  int neighbor_;
  PerturbRotate rotate_;
  bool inside_;

  // temporary
  PerturbTranslate translate_;
  Position displace_;
};

inline std::shared_ptr<PerturbMoveAVB> MakePerturbMoveAVB(
    argtype args = argtype()) {
  return std::make_shared<PerturbMoveAVB>(args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_PERTURB_MOVE_AVB_H_
