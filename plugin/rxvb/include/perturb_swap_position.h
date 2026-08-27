
#ifndef FEASST_CHAIN_PERTURB_SWAP_POSITION_H_
#define FEASST_CHAIN_PERTURB_SWAP_POSITION_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

class Position;
class RotationMatrix;

/**
 * Rigidly translate each of the two particles such that the positions of their
 * first sites are swapped, then randomly rotate each particle.
 */
class PerturbSwapPosition : public PerturbMove {
 public:
  explicit PerturbSwapPosition(argtype args = argtype());
  explicit PerturbSwapPosition(argtype * args);
  void move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbSwapPosition(std::istream& istr);
  virtual ~PerturbSwapPosition();

 protected:
  void serialize_perturb_position_swap_(std::ostream& ostr) const;

  // temporary and not serialized
  std::unique_ptr<Position> r0_, r1_, axis_, origin_, quat_;
  std::unique_ptr<RotationMatrix> rot_mat_;
};

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_SWAP_POSITION_H_
