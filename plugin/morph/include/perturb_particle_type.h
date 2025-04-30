
#ifndef FEASST_MORPH_PERTURB_PARTICLE_TYPE_H_
#define FEASST_MORPH_PERTURB_PARTICLE_TYPE_H_

#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

class Select;

/**
  Change the type of a Particle.
  If the Particle has multiple Sites, then fix the coordinate of the first site
  of both particle types to be equal,
  and then randomly rotate the rest of the sites according to the positions
  given in the fstprt file particle type definition (e.g., this only works
  properly for rigid molecules, not flexible ones).
 */
class PerturbParticleType : public Perturb {
 public:
  /**
    args:
    - type: type to set for particle.
   */
  explicit PerturbParticleType(argtype args = argtype());
  explicit PerturbParticleType(argtype * args);

  void set_particle_type(
    System * system,
    const Select& select,
    const int type);

  void perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held = false,
    Acceptance * acceptance = NULL) override;

  void revert(System * system) override;
  void finalize(System * system) override;
  std::string status_header() const override;
  std::string status() const override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbParticleType(std::istream& istr);
  virtual ~PerturbParticleType() {}

 private:
  int new_particle_type_;
  PerturbRotate rotate_;

  // temporary
  int old_particle_type_;
  Position tmp_pos_;
};

inline std::shared_ptr<PerturbParticleType> MakePerturbParticleType(
    argtype args = argtype()) {
  return std::make_shared<PerturbParticleType>(args);
}

}  // namespace feasst

#endif  // FEASST_MORPH_PERTURB_PARTICLE_TYPE_H_
