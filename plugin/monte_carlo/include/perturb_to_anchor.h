
#ifndef FEASST_MONTE_CARLO_PERTURB_TO_ANCHOR_H_
#define FEASST_MONTE_CARLO_PERTURB_TO_ANCHOR_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Assumes the selection is one site, and the anchor is one site.
  Move the selected site to the same position as the anchor.

  Note that this Perturb should always have 1 steps per stage,
  as its placement is completely deterministic.
 */
class PerturbToAnchor : public PerturbMove {
 public:
  explicit PerturbToAnchor(argtype args = argtype());
  explicit PerturbToAnchor(argtype * args);
  void move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbToAnchor(std::istream& istr);
  virtual ~PerturbToAnchor() {}

 protected:
  void serialize_perturb_to_anchor_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbToAnchor> MakePerturbToAnchor(
    argtype args = argtype()) {
  return std::make_shared<PerturbToAnchor>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TO_ANCHOR_H_
