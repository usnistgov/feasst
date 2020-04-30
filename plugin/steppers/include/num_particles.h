
#ifndef FEASST_STEPPERS_NUM_PARTICLES_H_
#define FEASST_STEPPERS_NUM_PARTICLES_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Accumulate statistics for number of particles.
 */
class NumParticles : public Analyze {
 public:
  /**
    args:
    - particle_type: index of particle type from configuration.
      If -1, sum all particles (default: -1).
    - group: index of group from configuration.
      If -1, sum all particles (default: -1).
      Can only be specified if particle_type is -1.
   */
  explicit NumParticles(const argtype &args = argtype());

  void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  const Accumulator& num_particles() const { return accumulator(); }

  // serialize
  std::string class_name() const override {
    return std::string("NumParticles"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<NumParticles>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit NumParticles(std::istream& istr);

 private:
  int particle_type_;
  int group_;
};

inline std::shared_ptr<NumParticles> MakeNumParticles(
    const argtype &args = argtype()) {
  return std::make_shared<NumParticles>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_NUM_PARTICLES_H_
