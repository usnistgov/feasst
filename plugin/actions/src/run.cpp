#include "utils/include/arguments.h"
#include "utils/include/arguments_extra.h"  // parse
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/system.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/modify.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "actions/include/run.h"

namespace feasst {

Run::Run(argtype * args) {
  class_name_ = "Run";
  num_trials_ = integer64("num_trials", args, -1);
  until_num_particles_ = integer("until_num_particles", args, -1);
  if (used("configuration_index", *args)) {
    WARN("Deprecated Run::configuration_index->config (see Configuration::name)");
  }
  configuration_index_ = integer("configuration_index", args, 0);
  config_ = str("config", args, "");
  particle_type_ = str("particle_type", args, "");
  for_hours_ = dble("for_hours", args, -1);
  until_criteria_complete_ = boolean("until_criteria_complete", args, false);
  if (str("until", args, "") == "complete") {
    until_criteria_complete_ = true;
  }
  until_file_exists_ = str("until_file_exists", args, "");
  trials_per_file_check_ = integer("trials_per_file_check", args, 1e5);
  const std::string trial_name_ = str("Trial", args, "");
  if (!trial_name_.empty()) {
    DEBUG("trial_name: " << trial_name_);
    if (!particle_type_.empty()) {
      args->insert({"particle_type", particle_type_});
    }
    if (!config_.empty()) {
      args->insert({"config", config_});
    }
    trial_ = parse(dynamic_cast<Trial*>(MakeTrial().get()), args, trial_name_, false);
    if (!trial_) {
      DEBUG("not a trial, try a factory");
      trials_ = parse(dynamic_cast<TrialFactoryNamed*>(std::make_shared<TrialFactoryNamed>().get()), args, trial_name_, false);
      ASSERT(trials_, "trial_name_:" << trial_name_ << " not found.");
    }
    DEBUG("extracting particle type");
    // if unused, extract particle_type
    // (e.g., TrialGrowFile gets particle_type from file)
    str("particle_type", args, "");
  }
  for (const std::string& stepper_name : split(str("Stepper", args, ""), ',')) {
    DEBUG("stepper_name: " << stepper_name);
    if (!config_.empty()) {
      args->insert({"config", config_});
    }
    std::shared_ptr<Analyze> an = parse(
      dynamic_cast<Analyze*>(std::make_shared<Analyze>().get()),
      args, stepper_name, false);
    if (an) {
      analyzers_.push_back(an);
    } else {
      std::shared_ptr<Modify> mod = parse(
        dynamic_cast<Modify*>(std::make_shared<Modify>().get()),
        args, stepper_name, false);
      if (mod) {
        modifiers_.push_back(mod);
      } else {
        FATAL("unrecognized stepper_name:" << stepper_name);
      }
    }
  }
}
Run::Run(argtype args) : Run(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(Run,);

Run::Run(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 3854 && version <= 3858, "mismatch version: " << version);
  feasst_deserialize(&num_trials_, istr);
  feasst_deserialize(&until_num_particles_, istr);
  if (version >= 3855) {
    feasst_deserialize(&configuration_index_, istr);
  }
  if (version >= 3857) {
    feasst_deserialize(&config_, istr);
  }
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&for_hours_, istr);
  feasst_deserialize(&until_criteria_complete_, istr);
  if (version >= 3856) {
    feasst_deserialize(&until_file_exists_, istr);
    feasst_deserialize(&trials_per_file_check_, istr);
  }
  if (version >= 3858) {
    feasst_deserialize(&num_trials_before_, istr);
    feasst_deserialize(&num_analyze_before_, istr);
    feasst_deserialize(&num_modify_before_, istr);
  }
}

void Run::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(3858, ostr);
  feasst_serialize(num_trials_, ostr);
  feasst_serialize(until_num_particles_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(config_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(for_hours_, ostr);
  feasst_serialize(until_criteria_complete_, ostr);
  feasst_serialize(until_file_exists_, ostr);
  feasst_serialize(trials_per_file_check_, ostr);
  feasst_serialize(num_trials_before_, ostr);
  feasst_serialize(num_analyze_before_, ostr);
  feasst_serialize(num_modify_before_, ostr);
}

void Run::run(MonteCarlo * mc) {
  if (!config_.empty()) {
    configuration_index_ = mc->system().configuration_index(config_);
  }
  if (trial_) {
    num_trials_before_ = mc->trials().num();
    mc->add(trial_);
  }
  if (trials_) {
    num_trials_before_ = mc->trials().num();
    mc->add(trials_);
  }
  if (analyzers_.size() > 0) {
    num_analyze_before_ = mc->num_analyzers();
    for (const std::shared_ptr<Analyze>& an : analyzers_) {
      mc->add(an);
    }
  }
  if (modifiers_.size() > 0) {
    num_modify_before_ = mc->num_modifiers();
    for (const std::shared_ptr<Modify>& mod : modifiers_) {
      mc->add(mod);
    }
  }
  mc->run_num_trials(num_trials_);
  mc->run_until_num_particles(until_num_particles_,
                              particle_type_,
                              configuration_index_);
  mc->run_for_hours(for_hours_);
  if (until_criteria_complete_) {
    mc->run_until_complete();
  }
  mc->run_until_file_exists(until_file_exists_, trials_per_file_check_);
  if (num_trials_before_ != -1) {
    for (int index = mc->trials().num() - 1; index >= num_trials_before_; --index) {
      mc->remove_trial(index);
    }
  }
  if (num_analyze_before_ != -1) {
    for (int index = mc->num_analyzers() - 1; index >= num_analyze_before_; --index) {
      mc->remove_analyze(index);
    }
  }
  if (num_modify_before_ != -1) {
    for (int index = mc->num_modifiers() - 1; index >= num_modify_before_; --index) {
      mc->remove_modify(index);
    }
  }
}

}  // namespace feasst
