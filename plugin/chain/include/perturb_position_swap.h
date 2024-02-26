
#ifndef FEASST_CHAIN_PERTURB_POSITION_SWAP_H_
#define FEASST_CHAIN_PERTURB_POSITION_SWAP_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Assumes the selection is two mobiles sites.
  Swap the positions of these.

  Note that this Perturb should always have 1 steps per stage,
  as its placement is completely deterministic.
 */
class PerturbPositionSwap : public PerturbMove {
 public:
  explicit PerturbPositionSwap(argtype args = argtype());
  explicit PerturbPositionSwap(argtype * args);
  void move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbPositionSwap(std::istream& istr);
  virtual ~PerturbPositionSwap() {}

 protected:
  void serialize_perturb_position_swap_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbPositionSwap> MakePerturbPositionSwap(
    argtype args = argtype()) {
  return std::make_shared<PerturbPositionSwap>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_POSITION_SWAP_H_
