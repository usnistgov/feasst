
#ifndef FEASST_GROWTH_EXPANDED_PERTURB_PARTICLE_TYPE_H_
#define FEASST_GROWTH_EXPANDED_PERTURB_PARTICLE_TYPE_H_

#include "monte_carlo/include/perturb.h"

namespace feasst {

/**
  Change the type of a site.
 */
class PerturbParticleType : public Perturb {
 public:
  /**
    args:
    - type: type to set for particle.
   */
  explicit PerturbParticleType(const argtype& args = argtype());

  void set_particle_type(
    System * system,
    const Select& select,
    const int type);

  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false) override;

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

  // temporary
  int old_particle_type_;
};

inline std::shared_ptr<PerturbParticleType> MakePerturbParticleType(const argtype& args = argtype()) {
  return std::make_shared<PerturbParticleType>(args);
}

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_PERTURB_PARTICLE_TYPE_H_
