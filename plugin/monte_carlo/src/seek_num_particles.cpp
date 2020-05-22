#include "utils/include/debug.h"
#include "utils/include/utils_io.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

SeekNumParticles::SeekNumParticles(const int num, const argtype& args) :
  num_(num) {
  Arguments args_(args);
  particle_type_ = args_.key("particle_type").dflt("0").integer();
  max_attempts_ = args_.key("max_attempts").dflt(str(1e8)).integer();
}

SeekNumParticles SeekNumParticles::with_metropolis(const argtype& args) {
  criteria_ = MakeMetropolis(args);
  return *this;
}

SeekNumParticles SeekNumParticles::with_trial_add(const argtype& args) {
  argtype typed_args;
  auto pair = args.find("particle_type");
  if (pair == args.end()) {
    typed_args.insert({"particle_type", str(particle_type_)});
  } else {
    ASSERT(pair->second == str(particle_type_), "particle type given to add: "
      << pair->second << " not equal to type of seek: " << particle_type_);
  }
  add_ = MakeTrialAdd(typed_args);
  return *this;
}

void SeekNumParticles::run(MonteCarlo * monte_carlo) {
  const Configuration& config = monte_carlo->configuration();
  ASSERT(config.num_particles_of_type(particle_type_) <= num_,
    "There are " << config.num_particles_of_type(particle_type_)
    << " particles of type " << particle_type_ << " which is > "
    << num_);
  Criteria * crit;
  if (criteria_) {
    crit = criteria_.get();
  } else {
    crit = monte_carlo->get_criteria();
  }
  System * sys = monte_carlo->get_system();
  Random * ran = monte_carlo->get_random();
  if (add_) {
    add_->precompute(crit, sys);
  }
  monte_carlo->reset_trial_stats();
  while (config.num_particles_of_type(particle_type_) < num_) {
    monte_carlo->get_trial_factory()->attempt(crit, sys, ran);
    if (add_) {
      add_->attempt(crit, sys, ran);
    }
    ASSERT(monte_carlo->trials().num_attempts() < max_attempts_,
      "max attempts:  "<< max_attempts_ << " reached during seek");
  }
  monte_carlo->initialize_energy();
  monte_carlo->reset_trial_stats();
}

}  // namespace feasst
