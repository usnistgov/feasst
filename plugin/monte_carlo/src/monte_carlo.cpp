#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

void MonteCarlo::add(const Configuration& config) {
  system_.add(config);
  config_set_ = true;
  if (potential_set_) system_set_ = true;
  ASSERT(!criteria_set_, "add config before criteria");
}

void MonteCarlo::add(const Potential& potential) {
  ASSERT(!criteria_set_, "add potential before criteria");
  system_.add(potential);
  system_.precompute();
  potential_set_ = true;
  if (config_set_) system_set_ = true;
}

void MonteCarlo::set(const int index, const Potential& potential) {
  // ASSERT(!criteria_set_, "add potential before criteria");
  ASSERT(potential_set_ || system_set_, "add potential before setting one");
  system_.set_unoptimized(index, potential);
  system_.precompute();
}

void MonteCarlo::set(const System& system) {
  system_set_ = true;
  system_ = system;
  system_.precompute();
  ASSERT(!criteria_set_, "add system before criteria");
}

void MonteCarlo::add(std::shared_ptr<Criteria> criteria) {
  ASSERT(system_set_, "set System before Criteria.");
  criteria->set_current_energy(system_.unoptimized_energy());
  DEBUG("current energy: " << criteria->current_energy());
  criteria_ = criteria;
  criteria_set_ = true;
}

void MonteCarlo::add(std::shared_ptr<Trial> trial) {
  ASSERT(criteria_set_, "set Criteria before Trials.");
  trial->precompute(criteria_.get(), &system_);
  trial_factory_.add(trial);
  // If later, perhaps after some initialization, more trials are added,
  // then Analyze and Modify classes may need to be re-initialized.
  // analyze_factory_.initialize(criteria_, system_, trial_factory_);
  // modify_factory_.initialize(criteria_, &system_, &trial_factory_);
}

void MonteCarlo::add(std::shared_ptr<Analyze> analyze) {
  ASSERT(criteria_set_, "set Criteria before Analyze");
  DEBUG("class name? " << analyze->class_name());
  if (analyze->is_multistate() and analyze->class_name() != "AnalyzeFactory") {
    auto multi = MakeAnalyzeFactory({{"multistate", "true"}});
    DEBUG("making multi " << multi->is_multistate());
    for (int state = 0; state < criteria_->num_states(); ++state) {
      DEBUG("state: " << state);
      std::shared_ptr<Analyze> an = deep_copy_derived(analyze);
      { std::stringstream ss;
        an->serialize(ss);
        DEBUG(ss.str());
      }
      an->set_state(state);
      // an->initialize(criteria_, system_, trial_factory_);
      multi->add(an);
    }
    analyze = multi;
  }
  analyze->initialize(criteria_.get(), &system_, &trial_factory_);
  DEBUG("mults " << analyze->is_multistate() << " class name? " << analyze->class_name());
  analyze_factory_.add(analyze);
}

void MonteCarlo::add(const std::shared_ptr<Modify> modify) {
  ASSERT(criteria_set_, "set Criteria before Modify");
  modify->initialize(criteria_.get(), &system_, &trial_factory_);
  modify_factory_.add(modify);
}

void MonteCarlo::seek_num_particles(const int num) {
  ASSERT(system_.configuration().num_particles() <= num,
    "assumes you only want to add particles, not delete them");
  auto add = MakeTrialAdd({{"particle_type", "0"}});
  add->precompute(criteria_.get(), &system_);
  while (system_.configuration().num_particles() < num) {
    attempt();
    add->attempt(criteria_.get(), &system_, random_.get());
  }
  trial_factory_.reset_stats();
}

}  // namespace feasst
