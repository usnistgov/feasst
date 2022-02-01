#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/progress_report.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

SeekNumParticles::SeekNumParticles(const int num, argtype args) :
  num_(num) {
  particle_type_ = integer("particle_type", &args, -1);
  max_attempts_ = integer("max_attempts", &args, 1e8);
  check_all_used(args);
}

SeekNumParticles SeekNumParticles::with_thermo_params(argtype args) {
  thermo_params_ = MakeThermoParams(args);
  return *this;
}

SeekNumParticles SeekNumParticles::with_metropolis() {
  criteria_ = MakeMetropolis();
  return *this;
}

SeekNumParticles SeekNumParticles::with_metropolis(
    std::shared_ptr<Constraint> constraint) {
  criteria_ = MakeMetropolis(constraint);
  return *this;
}

SeekNumParticles SeekNumParticles::with_trial_add(argtype args) {
  argtype typed_args = args;
  if (used("particle_type", args)) {
    const int ptype = integer("particle_type", &args);
    ASSERT(ptype == particle_type_, "particle type given to add: "
      << ptype << " not equal to type of seek: " << particle_type_);
  } else {
    std::string type = str(particle_type_);
    if (particle_type_ == -1) type = "0";
    typed_args.insert({"particle_type", type});
  }
  extra_trials_.add(MakeTrialAdd(typed_args));
  check_all_used(args);
  return *this;
}

SeekNumParticles SeekNumParticles::add(std::shared_ptr<Trial> trial) {
  extra_trials_.add(trial);
  return *this;
}

void SeekNumParticles::run(MonteCarlo * monte_carlo) {
  const Configuration& config = monte_carlo->configuration();
  ASSERT(config.num_particles_of_type(particle_type_) <= num_,
    "There are " << config.num_particles_of_type(particle_type_)
    << " particles of type " << particle_type_ << " which is > "
    << num_);
  std::shared_ptr<ThermoParams> original_params;
  if (thermo_params_) {
    if (monte_carlo->system().thermo_params_ptr_()) {
      original_params = std::make_shared<ThermoParams>(
        monte_carlo->system().thermo_params());
    }
    monte_carlo->set(thermo_params_);
  }
  Criteria * crit;
  if (criteria_) {
    crit = criteria_.get();
  } else {
    crit = monte_carlo->get_criteria();
  }
  ASSERT(crit, "SeekNumParticles requires criteria in MonteCarlo or input");
  System * sys = monte_carlo->get_system();
  crit->set_current_energy(sys->energy());
  Random * ran = monte_carlo->get_random();
  extra_trials_.precompute(crit, sys);
  monte_carlo->reset_trial_stats();
  int current_num = config.num_particles_of_type(particle_type_);
  if (report_) report_->set_num(num_ - current_num);
  while (current_num < num_) {
    int previous_num = current_num;
    DEBUG("current_num " << current_num);
    monte_carlo->get_trial_factory()->attempt(crit, sys, ran);
    extra_trials_.attempt(crit, sys, ran);
    ASSERT(monte_carlo->trials().num_attempts() < max_attempts_,
      "max attempts:  "<< max_attempts_ << " reached during seek");
    current_num = config.num_particles_of_type(particle_type_);
    if (report_) {
      for (int diff = 0; diff < current_num - previous_num; ++diff) {
        report_->check();
      }
    }
  }
  if (thermo_params_ && original_params) monte_carlo->set(original_params);
  monte_carlo->initialize_criteria();
  monte_carlo->reset_trial_stats();
}

}  // namespace feasst
