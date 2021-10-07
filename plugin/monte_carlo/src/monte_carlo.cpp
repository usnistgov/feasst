#include <iostream>
#include <fstream>
#include "utils/include/serialize.h"
#include "utils/include/checkpoint.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/action.h"

// for parsing factories
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

MonteCarlo::MonteCarlo(std::shared_ptr<Random> random) {
  set(random);
//    timer_other_ = timer_.add("other");
//    timer_trial_ = timer_.add("trial");
//    timer_analyze_ = timer_.add("analyze");
//    timer_modify_ = timer_.add("modify");
//    timer_checkpoint_ = timer_.add("checkpoint");
}

MonteCarlo::MonteCarlo() : MonteCarlo(std::make_shared<RandomMT19937>()) {}

void MonteCarlo::parse_(arglist * args) {
  DEBUG("first " << args->begin()->first);
  std::cout << "{" << args->begin()->first << ","
            << str(args->begin()->second) << "}," << std::endl;

  // parse all derived classes of Random
  std::shared_ptr<Random> ran =
    parse(dynamic_cast<Random*>(MakeRandomMT19937().get()), args);
  if (ran) {
    DEBUG("parsing Random");
    set(ran);
    return;
  }

  // parse Checkpoint
  if (args->begin()->first == "Checkpoint") {
    DEBUG("parsing Checkpoint");
    set(MakeCheckpoint(args->begin()->second));
    args->erase(args->begin());
    return;
  }

  // parse Configuration
  if (args->begin()->first == "Configuration") {
    DEBUG("parsing Configuration");
    add(MakeConfiguration(args->begin()->second));
    args->erase(args->begin());
    return;
  }

  // parse Potential
  if (args->begin()->first == "Potential") {
    DEBUG("parsing Potential");
    add(MakePotential(args->begin()->second));
    args->erase(args->begin());
    return;
  }

  // parse ThermoParams
  if (args->begin()->first == "ThermoParams") {
    DEBUG("parsing ThermoParams");
    set(MakeThermoParams(args->begin()->second));
    args->erase(args->begin());
    return;
  }

  // parse all derived classes of Criteria
  std::shared_ptr<Criteria> crit =
    parse(dynamic_cast<Criteria*>(MakeMetropolis().get()), args);
  if (crit) {
    DEBUG("parsing Criteria");
    set(crit);
    return;
  }

  // parse all derived classes of Trial
  std::shared_ptr<Trial> trial =
    parse(dynamic_cast<Trial*>(MakeTrial().get()), args);
  if (trial) {
    DEBUG("parsing Trial");
    add(trial);
    return;
  }

  // parse all derived classes of Analyze
  std::shared_ptr<Analyze> an =
    parse(dynamic_cast<Analyze*>(std::make_shared<Analyze>().get()), args);
  if (an) {
    DEBUG("parsing Analyze");
    add(an);
    return;
  }

  // parse all derived classes of Modify
  std::shared_ptr<Modify> mod =
    parse(dynamic_cast<Modify*>(std::make_shared<Modify>().get()), args);
  if (mod) {
    DEBUG("parsing Modify");
    add(mod);
    return;
  }

  // parse all derived classes of Action
  std::shared_ptr<Action> act =
    parse(dynamic_cast<Action*>(std::make_shared<Action>().get()), args);
  if (act) {
    DEBUG("parsing Action");
    run(act);
    return;
  }

}

MonteCarlo::MonteCarlo(arglist args) : MonteCarlo() {
  args_ = args;
  int size = static_cast<int>(args_.size());
  int previous_size = size;
  while (size > 0) {
    previous_size = size;
    DEBUG("size " << size);
    parse_(&args_);
    size = static_cast<int>(args_.size());
    ASSERT(previous_size - 1 == size,
      "Unrecognized argument: " << args_.begin()->first);
  }
}

void MonteCarlo::run(std::shared_ptr<Action> action) {
  action_ = action;
  action_->run(this);
//  action_->run(&system_, criteria_, &trial_factory_, &analyze_factory_,
//               &modify_factory_, checkpoint_, random_);
}

void MonteCarlo::seed_random(const int seed) {
  random_->seed(seed);
}

void MonteCarlo::add(std::shared_ptr<Configuration> config) {
  system_.add(*config);
  config_set_ = true;
  if (potential_set_) system_set_ = true;
  ASSERT(!criteria_set_, "add config before criteria");
}

void MonteCarlo::add(const Configuration& config) {
  WARN("Use MakeConfiguration instead of Configuration");
  system_.add(config);
  config_set_ = true;
  if (potential_set_) system_set_ = true;
  ASSERT(!criteria_set_, "add config before criteria");
}

void MonteCarlo::add(std::shared_ptr<Potential> potential) {
  ASSERT(!criteria_set_, "add potential before criteria");
  ASSERT(config_set_ || system_set_, "config:" << config_set_ <<
    " or system:" << system_set_ << " must be set before adding a potential");
  system_.add(potential);
  system_.precompute();
  potential_set_ = true;
}

void MonteCarlo::set(const int index, std::shared_ptr<Potential> potential) {
  // ASSERT(!criteria_set_, "add potential before criteria");
  ASSERT(potential_set_ || system_set_, "add potential before setting one");
  system_.set_unoptimized(index, potential);
  system_.precompute();
}

void MonteCarlo::set(std::shared_ptr<ThermoParams> thermo_params) {
  system_.set(thermo_params);
  thermo_params_set_ = true;
  if (config_set_ && potential_set_) system_set_ = true;
}

void MonteCarlo::set(const System& system) {
  system_set_ = true;
  system_ = system;
  system_.precompute();
  // ASSERT(!criteria_set_, "add system before criteria");
  // HWH used in clones.cpp to transfer configurations
}

void MonteCarlo::set(std::shared_ptr<Criteria> criteria) {
  ASSERT(system_set_, "set System before Criteria.");
  criteria_ = criteria;
  criteria_set_ = true;
  initialize_criteria();
  // criteria->set_current_energy(system_.unoptimized_energy());
  DEBUG("current energy: " << criteria->current_energy());
}

void MonteCarlo::add(std::shared_ptr<Trial> trial) {
  ASSERT(criteria_set_, "set Criteria before Trials.");

  // Error check Ewald
  for (const std::shared_ptr<Potential> pot : system_.potentials().potentials()) {
    if (pot->visit_model().class_name() == "Ewald" ||
        pot->visit_model().class_name() == "LongRangeCorrections") {
      for (int stage = 0; stage < trial->num_stages(); ++stage) {
        // Require reference potentials for multi-stage trials.
        if (trial->num_stages() > 1 && trial->stage(stage).reference() == -1) {
          ERROR(trial->class_name() << " " << trial->description()
            << " should use a reference potential "
            << "without Ewald or LongRangeCorrections due to multiple stages "
            << "which complicate revert");
        }
      }
    }
  }

  // Error check DCCB
  for (int stage = 0; stage < trial->num_stages(); ++stage) {
    if (trial->stage(stage).rosenbluth().num() > 1) {
      if (trial->description() == "TrialTranslate" ||
          trial->description() == "TrialRotate") {
        ERROR(trial->description() << " cannot be used with multiple steps.");
      }
    }
  }

  // Don't allow new_only with TrialTransfer.
  if (trial->num_stages() > 0) {
    if (trial->stage(0).is_new_only()) {
      if (trial->stage(0).perturb().class_name() == "PerturbRemove") {
        FATAL("PerturbRemove not implemented with new_only due to trial_state");
      }
    }
  }

  // flatten TrialFactory by adding the individual trials instead.
  if (trial->class_name() == "TrialFactory") {
    DEBUG("flattening");
    double total_weight = 0.;
    for (std::shared_ptr<Trial> itrl : trial->trials()) {
      total_weight += itrl->weight();
    }
    for (std::shared_ptr<Trial> itrl : trial->trials()) {
      itrl->set_weight(trial->weight()*itrl->weight()/total_weight);
      DEBUG(itrl->weight() << " " << trial->weight() << " " << total_weight);
      add(itrl);
    }
  } else {
    trial->precompute(criteria_.get(), &system_);
    trial_factory_.add(trial);
  }

  // HWH depreciate?
  // If later, perhaps after some initialization, more trials are added,
  // then Analyze and Modify classes may need to be re-initialized.
  // analyze_factory_.initialize(criteria_, system_, trial_factory_);
  // modify_factory_.initialize(criteria_, &system_, &trial_factory_);
}

bool MonteCarlo::duplicate_stepper_file_name_(const std::string file_name) {
  if (!file_name.empty()) {
    for (const std::shared_ptr<Analyze> an : analyze_factory_.analyzers()) {
      if (an->file_name() == file_name) return true;
    }
    for (const std::shared_ptr<Modify> mod : modify_factory_.modifiers()) {
      if (mod->file_name() == file_name) return true;
    }
  }
  return false;
}

void MonteCarlo::add(std::shared_ptr<Analyze> analyze) {
  ASSERT(criteria_set_, "set Criteria before Analyze");
  ASSERT(!duplicate_stepper_file_name_(analyze->file_name()),
    "Analyze " << analyze->class_name() << " should not have the same " <<
    "file_name as an already existing Stepper file name");

  // process multistate
  DEBUG("class name? " << analyze->class_name());
  if (analyze->is_multistate() && analyze->class_name() != "AnalyzeFactory") {
    int steps_per_write = 1;
    std::string file_name;
    if (analyze->is_multistate_aggregate()) {
      steps_per_write = analyze->steps_per_write();
      file_name = analyze->file_name();
      analyze->initialize(criteria_.get(), &system_, &trial_factory_);
    }
    auto multi = MakeAnalyzeFactory({
      {"multistate", "true"},
      {"steps_per_write", str(steps_per_write)},
      {"file_name", file_name},
      {"steps_per_update", "1"},
      {"stop_after_phase", str(analyze->stop_after_phase())},
      {"start_after_phase", str(analyze->start_after_phase())},
      {"file_name_append_phase", str(analyze->file_name_append_phase())},
      {"multistate_aggregate", str(analyze->is_multistate_aggregate())}});
    DEBUG("making multi " << multi->is_multistate());
    for (int state = 0; state < criteria_->num_states(); ++state) {
      DEBUG("state: " << state);
      std::shared_ptr<Analyze> an = deep_copy_derived(analyze);
      { std::stringstream ss;
        an->serialize(ss);
        DEBUG(ss.str());
      }
      if (analyze->is_multistate_aggregate()) {
        an->empty_file_name();
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

// copied from above
void MonteCarlo::add(std::shared_ptr<Modify> modify) {
  ASSERT(criteria_set_, "set Criteria before Modify");
  ASSERT(!duplicate_stepper_file_name_(modify->file_name()),
    "Modify " << modify->class_name() << " should not have the same " <<
    "file_name as an already existing Stepper file name");
  DEBUG("class name? " << modify->class_name());
  if (modify->is_multistate() && modify->class_name() != "ModifyFactory") {
    int steps_per_write = 1;
    std::string file_name;
    DEBUG("aggregate? " << modify->is_multistate_aggregate());
    if (modify->is_multistate_aggregate()) {
      steps_per_write = modify->steps_per_write();
      file_name = modify->file_name();
      modify->initialize(criteria_.get(), &system_, &trial_factory_);
    }
    auto multi = MakeModifyFactory({
      {"multistate", "true"},
      {"steps_per_write", str(steps_per_write)},
      {"file_name", file_name},
      {"steps_per_update", "1"},
      {"stop_after_phase", str(modify->stop_after_phase())},
      {"start_after_phase", str(modify->start_after_phase())},
      {"file_name_append_phase", str(modify->file_name_append_phase())},
      {"multistate_aggregate", str(modify->is_multistate_aggregate())}});
    DEBUG("making multi " << multi->is_multistate());
    for (int state = 0; state < criteria_->num_states(); ++state) {
      DEBUG("state: " << state);
      std::shared_ptr<Modify> mod = deep_copy_derived(modify);
      { std::stringstream ss;
        mod->serialize(ss);
        DEBUG(ss.str());
      }
      if (modify->is_multistate_aggregate()) {
        mod->empty_file_name();
      }
      mod->set_state(state);
      // mod->initialize(criteria_.get(), &system_, &trial_factory_);
      multi->add(mod);
    }
    modify = multi;
  }
  modify->initialize(criteria_.get(), &system_, &trial_factory_);
  DEBUG("mults " << modify->is_multistate() << " class name? " << modify->class_name());
  modify_factory_.add(modify);
}

//void MonteCarlo::add(const std::shared_ptr<Modify> modify) {
//  ASSERT(criteria_set_, "set Criteria before Modify");
//  modify->initialize(criteria_.get(), &system_, &trial_factory_);
//  modify_factory_.add(modify);
//}

void MonteCarlo::set(const std::shared_ptr<Checkpoint> checkpoint) {
  checkpoint_ = checkpoint;
}

void MonteCarlo::after_trial_modify_() {
  modify_factory_.trial(criteria_.get(), &system_, &trial_factory_);
  if (checkpoint_) {
    checkpoint_->check(*this);
  }
}

void MonteCarlo::serialize(std::ostream& ostr) const {
  feasst_serialize_version(529, ostr);
  feasst_serialize_fstobj(system_, ostr);
  feasst_serialize_fstdr(criteria_, ostr);
  feasst_serialize_fstobj(trial_factory_, ostr);
  feasst_serialize_fstobj(analyze_factory_, ostr);
  feasst_serialize_fstobj(modify_factory_, ostr);
  feasst_serialize(checkpoint_, ostr);
  feasst_serialize_fstdr(random_, ostr);
  feasst_serialize_fstdr(action_, ostr);
  feasst_serialize(args_, ostr);
  feasst_serialize(config_set_, ostr);
  feasst_serialize(potential_set_, ostr);
  feasst_serialize(thermo_params_set_, ostr);
  feasst_serialize(system_set_, ostr);
  feasst_serialize(criteria_set_, ostr);
  feasst_serialize_endcap("MonteCarlo", ostr);
}

MonteCarlo::MonteCarlo(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 529, "version: " << version);
  feasst_deserialize_fstobj(&system_, istr);
  // feasst_deserialize_fstdr(criteria_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      criteria_ = criteria_->deserialize(istr);
    }
  }
  feasst_deserialize_fstobj(&trial_factory_, istr);
  feasst_deserialize_fstobj(&analyze_factory_, istr);
  feasst_deserialize_fstobj(&modify_factory_, istr);
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(checkpoint_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      checkpoint_ = std::make_shared<Checkpoint>(istr);
    }
  }
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize_fstdr(random_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      random_ = random_->deserialize(istr);
    }
  }
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize_fstdr(action_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      action_ = action_->deserialize(istr);
    }
  }
  feasst_deserialize(&args_, istr);
  feasst_deserialize(&config_set_, istr);
  feasst_deserialize(&potential_set_, istr);
  feasst_deserialize(&thermo_params_set_, istr);
  feasst_deserialize(&system_set_, istr);
  feasst_deserialize(&criteria_set_, istr);
  feasst_deserialize_endcap("MonteCarlo", istr);
}

void MonteCarlo::load_cache_(const bool load) {
  random_->set_cache_to_load(load);
  system_.load_cache(load);
}

void MonteCarlo::unload_cache_(const MonteCarlo& mc) {
  random_->set_cache_to_unload((*mc.random_));
  system_.unload_cache(mc.system());
}

void MonteCarlo::before_attempts_() {
  ASSERT(system_set_, "system must be set before attempting trials.");
  ASSERT(criteria_set_, "criteria must be set before attempting trials.");
}

void MonteCarlo::revert_(const int trial_index,
    const bool accepted,
    const bool auto_reject,
    const double ln_prob) {
  trial_factory_.revert(trial_index, accepted, auto_reject, &system_);
  DEBUG("reverting " << criteria_->current_energy());
  criteria_->revert_(accepted, ln_prob);
}

void MonteCarlo::attempt_(int num_trials,
    TrialFactory * trial_factory,
    Random * random) {
  before_attempts_();
  for (int trial = 0; trial < num_trials; ++trial) {
    DEBUG("mc trial: " << trial);
    trial_factory->attempt(criteria_.get(), &system_, random);
    after_trial_analyze_();
    after_trial_modify_();
  }
}

bool MonteCarlo::attempt_trial(const int index) {
  return trial_factory_.attempt(criteria_.get(), &system_,
                                index, random_.get());
}

void MonteCarlo::imitate_trial_rejection_(const int trial_index,
    const double ln_prob,
    const bool auto_reject,
    const int state_old,
    const int state_new) {
  trial_factory_.imitate_trial_rejection_(trial_index, auto_reject);
  criteria_->imitate_trial_rejection_(ln_prob, state_old, state_new);
}

double MonteCarlo::initialize_system() {
  system_.precompute();
  const double en = system_.unoptimized_energy();
  system_.energy();
  for (int ref = 0; ref < system_.num_references(); ++ref) {
    system_.reference_energy(ref);
  }
  return en;
}

void MonteCarlo::initialize_criteria() {
  const double en = initialize_system();
  if (criteria_) {
    criteria_->set_current_energy(en);
  }
}

void MonteCarlo::initialize_trials() {
  for (int trial = 0; trial < trial_factory_.num(); ++trial) {
    trial_factory_.get_trial(trial)->precompute(criteria_.get(), &system_);
  }
}

void MonteCarlo::initialize_analyzers() {
  for (int an = 0; an < analyze_factory_.num(); ++an) {
    analyze_factory_.get_analyze(an)->initialize(
      criteria_.get(), &system_, &trial_factory_);
  }
}

void MonteCarlo::write_checkpoint() const {
  if (checkpoint_) checkpoint_->write(*this);
}

void MonteCarlo::run_until_complete_(TrialFactory * trial_factory,
                                     Random * random) {
  while (!criteria_->is_complete()) {
    attempt_(1, trial_factory, random);
  }
  write_checkpoint();
}

void MonteCarlo::synchronize_(const MonteCarlo& mc, const Select& perturbed) {
  system_.synchronize_(mc.system(), perturbed);
  criteria_->synchronize_(mc.criteria());
  trial_factory_.synchronize_(mc.trials());
}

std::shared_ptr<MonteCarlo> MakeMonteCarlo(const std::string file_name) {
  std::ifstream file(file_name);
  std::string line;
  std::getline(file, line);
  ASSERT(!line.empty(), "file: " << file_name << " is empty");
  std::stringstream ss;
  ss << line;
  return std::make_shared<MonteCarlo>(ss);
}

void MonteCarlo::add(const Potential& potential) {
  FATAL("Please use MakePotential instead of Potential");
}

}  // namespace feasst
