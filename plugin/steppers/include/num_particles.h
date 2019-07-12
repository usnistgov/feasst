
#ifndef FEASST_STEPPERS_NUM_PARTICLES_H_
#define FEASST_STEPPERS_NUM_PARTICLES_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate statistics for number of particles.
 */
class NumParticles : public Analyze {
 public:
  NumParticles(
    /**
      num_block : number of updated per block (default: 1e5).

      particle_type : index of particle type from configuration.
        If -1, sum all particles (default: -1).

      group : index of group from configuration.
        If -1, sum all particles (default: -1).
        Can only be specified if particle_type is -1.
     */
    const argtype &args = argtype());

  void update(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    ASSERT(particle_type_ == -1 || group_ == -1,
      "both particle type(" << particle_type_ << ") and group(" << group_ <<
      ") cannot be specified at the same time.");
    const Configuration& config = system.configuration();
    DEBUG(particle_type_);
    DEBUG(group_);
    if (particle_type_ == -1) {
      if (group_ == -1) {
        num_particles_.accumulate(config.num_particles(0));
      } else {
        num_particles_.accumulate(config.num_particles(group_));
      }
    } else {
      num_particles_.accumulate(config.num_particles_of_type(particle_type_));
    }
  }

  std::string write(const Criteria * criteria,
      const System& system,
      const TrialFactory& trial_factory) override {
    std::stringstream ss;
    ss << num_particles_.str() << " ";
    DEBUG(ss.str());
    return ss.str();
  }

  const Accumulator& num_particles() const { return num_particles_; }

  const Accumulator& accumulator() const override { return num_particles(); }

  std::string class_name() const override { return std::string("NumParticles"); }

  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<NumParticles>(istr); }

  void serialize(std::ostream& ostr) const override {
    Stepper::serialize(ostr);
    feasst_serialize_version(1, ostr);
    feasst_serialize_fstobj(num_particles_, ostr);
    feasst_serialize(particle_type_, ostr);
    feasst_serialize(group_, ostr);
  }

  NumParticles(std::istream& istr) : Analyze(istr) {
    feasst_deserialize_version(istr);
    feasst_deserialize_fstobj(&num_particles_, istr);
    feasst_deserialize(&particle_type_, istr);
    feasst_deserialize(&group_, istr);
  }

 private:
  Accumulator num_particles_;
  int particle_type_;
  int group_;
};

inline std::shared_ptr<NumParticles> MakeNumParticles(const argtype &args = argtype()) {
  return std::make_shared<NumParticles>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_NUM_PARTICLES_H_
