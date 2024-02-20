
#ifndef FEASST_CHAIN_PERTURB_PARTICLE_PIVOT_H_
#define FEASST_CHAIN_PERTURB_PARTICLE_PIVOT_H_

#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/**
  Rotate the positions of selection about the anchor site.
 */
class PerturbParticlePivot : public PerturbRotate {
 public:
  explicit PerturbParticlePivot(argtype args = argtype());
  explicit PerturbParticlePivot(argtype * args);

  void move(const bool is_position_held, System * system, TrialSelect * select,
            Random * random, Acceptance * acceptance) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbParticlePivot(std::istream& istr);
  virtual ~PerturbParticlePivot() {}

 protected:
  void serialize_perturb_particle_pivot_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbParticlePivot> MakePerturbParticlePivot(
    argtype args = argtype()) {
  return std::make_shared<PerturbParticlePivot>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_PARTICLE_PIVOT_H_
