#include <iostream>
#include <fstream>
#include "utils/include/serialize_extra.h"
#include "utils/include/arguments.h"
#include "utils/include/file.h"
#include "utils/include/io.h"
#include "utils/include/checkpoint.h"
#include "utils/include/timer_rdtsc.h"
#include "configuration/include/neighbor_criteria.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/visit_model.h"
#include "system/include/potential.h"
#include "system/include/system.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/rosenbluth.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/action.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/monte_carlo.h"

// for parsing factories
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

/// If args contains derived class of T, return factory pointer and remove from
/// args.
template <class T>
std::shared_ptr<T> parse(T * obj, arglist * args) {
  std::shared_ptr<T> new_obj;
  const auto& map = obj->deserialize_map();
  DEBUG("parsing " << args->begin()->first);
  if (map.count(args->begin()->first) > 0) {
    new_obj = obj->factory(args->begin()->first, &args->begin()->second);
    DEBUG(new_obj->class_name());
    feasst_check_all_used(args->begin()->second);
    return new_obj;
  }
  return new_obj;
}

void MonteCarlo::set_timer() {
  timer_ = std::make_unique<TimerRDTSC>(5);
  timer_->start(4);
}

MonteCarlo::MonteCarlo(std::shared_ptr<Random> random) {
  set(random);
  system_ = std::make_unique<System>();
  trial_factory_ = std::make_unique<TrialFactory>();
  analyze_factory_ = std::make_unique<AnalyzeFactory>();
  modify_factory_ = std::make_unique<ModifyFactory>();
}

MonteCarlo::MonteCarlo() : MonteCarlo(std::make_shared<RandomMT19937>()) {}
MonteCarlo::~MonteCarlo() {}

void MonteCarlo:: record_next_arg_(arglist *args) {
  if (args->size() > 0) {
    next_arg_ = *args->begin();
  }
}

void MonteCarlo::parse_args(arglist * args, const bool silent) {
  DEBUG("first " << args->begin()->first);

  // Check for deprecated names
  if (args->begin()->first == "PressureFromTestVolume") {
    WARN("PressureFromTestVolume was renamed to GhostTrialVolume");
    args->begin()->first = "GhostTrialVolume";
  } else if (args->begin()->first == "ProfileTrials") {
    WARN("ProfileTrials is deprecated. Use ProfileCPU instead.");
  } else if (args->begin()->first == "RemoveTrial") {
    WARN("RemoveTrial is deprecated. Use Remove instead.");
  } else if (args->begin()->first == "RemoveAnalyze") {
    WARN("RemoveAnalyze is deprecated. Use Remove instead.");
  } else if (args->begin()->first == "RemoveModify") {
    WARN("RemoveModify is deprecated. Use Remove instead.");
  }

  // repeat for each copy (which may be a config)
  DEBUG("parse_for_num_configs " << parse_for_num_configs_);
  ASSERT(parse_for_num_configs_ == 1 || parse_for_num_configs_ == 2,
    "hard corded for 1 or two configs");
  int num_copy = parse_for_num_configs_;
  DEBUG("num_copy " << num_copy);
  for (int copy = 0; copy <= num_copy; ++copy) {
    DEBUG("copy " << copy);
    if (parse_for_num_configs_ > 1) {
      if (copy == 0 && args->begin()->first != "EndCopy") {
        args->begin()->second.insert({"configuration_index", str(copy)});
        args->insert(args->begin() + 1, *args->begin());
      } else if (copy == 1 && args->begin()->first != "EndCopy") {
        args->begin()->second.insert({"configuration_index", str(copy)});
        auto pair = args->begin()->second.find("configuration_index");
        ASSERT(pair != args->begin()->second.end(), "err");
        pair->second = str(copy);
        //for (argtype::iterator it = args->begin()->second.begin();
        //     it != args->begin()->second.end(); ++it ) {
        for (auto &p : args->begin()->second) {
          for (const std::vector<std::string>& vals : parse_replace_) {
            feasst::replace(vals[0], vals[1], &p.second);
          }
        }
      }
      for (auto &p : args->begin()->second) {
        if (!replace_with_index_.empty()) {
          feasst::replace(replace_with_index_, str(copy), &p.second);
        }
      }
    }
    if (!silent &&
        args->begin()->first != "CopyFollowingLines" &&
        args->begin()->first != "EndCopy") {
      std::cout << args->begin()->first;
      std::string second = str(args->begin()->second);
      if (!second.empty() && second != "=") {
        std::cout << " " << second;
      }
      std::cout << std::endl;
    }

    // parse all derived classes of Random
    std::shared_ptr<Random> ran =
      parse(dynamic_cast<Random*>(MakeRandomMT19937().get()), args);
    if (ran) {
      DEBUG("parsing Random");
      set(ran);
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse Checkpoint
    if (args->begin()->first == "Checkpoint") {
      DEBUG("parsing Checkpoint");
      set(MakeCheckpoint(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse Configuration
    if (args->begin()->first == "Configuration") {
      DEBUG("parsing Configuration");
      add(MakeConfiguration(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse NeighborCriteria
    if (args->begin()->first == "NeighborCriteria") {
      DEBUG("parsing NeighborCriteria");
      add(MakeNeighborCriteria(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse Potential
    if (args->begin()->first == "Potential") {
      DEBUG("parsing Potential");
      //const int pconfig = integer("configuration_index", &(args->begin()->second), 0);
      add(MakePotential(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse reference Potential
    if (args->begin()->first == "ReferencePotential") {
      WARN("Deprecated ReferencePotential->RefPotential");
      DEBUG("parsing ReferencePotential");
      add_to_reference(MakePotential(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse optimized Potential
    if (args->begin()->first == "OptimizedPotential") {
      DEBUG("parsing OptimizedPotential");
      add_to_optimized(MakePotential(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse ThermoParams
    if (args->begin()->first == "ThermoParams") {
      DEBUG("parsing ThermoParams");
      set(MakeThermoParams(args->begin()->second));
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse all derived classes of Criteria
    std::shared_ptr<Criteria> crit =
      parse(dynamic_cast<Criteria*>(MakeMetropolis().get()), args);
    if (crit) {
      DEBUG("parsing Criteria");
      set(crit);
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse all derived classes of Trial
    std::shared_ptr<Trial> trial =
      parse(dynamic_cast<Trial*>(MakeTrial().get()), args);
    if (trial) {
      DEBUG("parsing Trial");
      add(trial);
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse all derived classes of TrialFactoryNamed
    std::shared_ptr<TrialFactoryNamed> trials =
      parse(dynamic_cast<TrialFactoryNamed*>(std::make_shared<TrialFactoryNamed>().get()), args);
    if (trials) {
      DEBUG("parsing TrialFactoryNamed");
      add(trials);
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse all derived classes of Analyze
    std::shared_ptr<Analyze> an =
      parse(dynamic_cast<Analyze*>(std::make_shared<Analyze>().get()), args);
    if (an) {
      DEBUG("parsing Analyze");
      add(an);
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse all derived classes of Modify
    std::shared_ptr<Modify> mod =
      parse(dynamic_cast<Modify*>(std::make_shared<Modify>().get()), args);
    if (mod) {
      DEBUG("parsing Modify");
      add(mod);
      args->erase(args->begin());
      if (copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }

    // parse all derived classes of Action
    std::shared_ptr<Action> act =
      parse(dynamic_cast<Action*>(std::make_shared<Action>().get()), args);
    if (act) {
      DEBUG("parsing Action");
      std::string aname = args->begin()->first;
      DEBUG("aname " << aname);
      args->erase(args->begin());
      record_next_arg_(args);
      run(act);
      if (aname == "CopyFollowingLines" ||
          aname == "EndCopy" || copy + 1 == num_copy) {
        return;
      } else {
        continue;
      }
    }
  }
}

void MonteCarlo::begin(arglist args, const bool silent) {
  args_ = args;
  resume(silent);
}

void MonteCarlo::resume(const bool silent) {
  if (action_) {
    run(action_);
  }
  int size = static_cast<int>(args_.size());
  DEBUG("size " << size);
  int previous_size = size;
  while (size > 0) {
    previous_size = size;
    DEBUG("size " << size);
    parse_args(&args_, silent);
    size = static_cast<int>(args_.size());
    ASSERT(previous_size - 1 == size,
      "Unrecognized argument: " << args_.begin()->first);
  }
}

MonteCarlo::MonteCarlo(arglist args, const bool silent) : MonteCarlo() {
  begin(args, silent);
}

void MonteCarlo::add_args(arglist args) {
  args_.insert(args_.end(), args.begin(), args.end());
}

void MonteCarlo::run(std::shared_ptr<Action> action) {
  action_ = action;
  action_->run(this);
//  action_->run(&system_, criteria_, trial_factory_.get(), analyze_factory_.get(),
//               modify_factory_.get(), checkpoint_, random_);
}

void MonteCarlo::seed_random(const int seed) {
  random_->seed(seed);
}

void MonteCarlo::add(std::shared_ptr<Configuration> config) {
  if (config->num_particle_types() == 0) {
    FATAL("There are no particle types in config");
  }
  system_->add(config);
  config_set_ = true;
  if (potential_set_) system_set_ = true;
  ASSERT(!criteria_set_, "add config before criteria");
}

void MonteCarlo::potential_check_(const Potential& pot) {
  // Check that if the Domain is tilted, then the cell list is only used in an
  // OptPotential to catch possible errors.
  if (pot.visit_model().class_name() == "VisitModelCell") {
    if (system_->configuration(pot.configuration_index()).domain().is_tilted()) {
      WARN("OptPotential with VisitModelCell is recommended to catch possible "
        << " errors if the Configuration/Domain is triclinic");
    }
  }
}

void MonteCarlo::add(std::shared_ptr<Potential> potential) {
  ASSERT(!criteria_set_, "add potential before criteria");
  ASSERT(config_set_ || system_set_, "config:" << config_set_ <<
    " or system:" << system_set_ << " must be set before adding a potential");
  system_->add(potential);
  potential_check_(*potential.get());
  system_->precompute();
  potential_set_ = true;
}

void MonteCarlo::set(const int index, std::shared_ptr<Potential> potential) {
  // ASSERT(!criteria_set_, "add potential before criteria");
  ASSERT(potential_set_ || system_set_, "add potential before setting one");
  ASSERT(system_->num_configurations() <= 1, "not implemented");
  system_->set_unoptimized(index, potential);
  potential_check_(*potential.get());
  system_->precompute();
}

void MonteCarlo::set(std::shared_ptr<ThermoParams> thermo_params) {
  system_->set(thermo_params);
  thermo_params_set_ = true;
  if (config_set_ && potential_set_) system_set_ = true;
}

void MonteCarlo::set(const System& system) {
  system_set_ = true;
  system_ = std::make_unique<System>(system);
  system_->precompute();
  // ASSERT(!criteria_set_, "add system before criteria");
  // HWH used in clones.cpp to transfer configurations
}

void MonteCarlo::set(std::shared_ptr<Criteria> criteria) {
  ASSERT(system_set_, "set System before Criteria.");
  criteria_ = criteria;
  criteria_set_ = true;
  initialize_criteria();
  // criteria->set_current_energy(system_->unoptimized_energy());
  DEBUG("current energy: " << criteria->current_energy());
}

void MonteCarlo::add(std::shared_ptr<Trial> trial) {
  ASSERT(criteria_set_, "set Criteria before Trials.");

  // Error check Ewald
  for (const std::shared_ptr<Potential>& pot : system_->potentials().potentials()) {
    if (pot->visit_model().class_name() == "Ewald" ||
        pot->visit_model().class_name() == "LongRangeCorrections") {
      for (int stage = 0; stage < trial->num_stages(); ++stage) {
        // Require reference potentials for multi-stage trials.
        const TrialStage& ts = trial->stage(stage);
        if (trial->num_stages() > 1 && (ts.reference() == -1 && ts.ref().empty())) {
          ERROR(trial->class_name() << " " << trial->description()
            << " should use a reference potential "
            << "without Ewald or LongRangeCorrections due to multiple stages "
            << "which complicate revert.");
        }
      }
      // Require reference potentials for any trial with multiple steps
      for (int stage = 0; stage < trial->num_stages(); ++stage) {
        const TrialStage& ts = trial->stage(stage);
        ASSERT(ts.num_steps() == 1 || ts.reference() != -1 || !ts.ref().empty(),
          "Ewald and LongRangeCorrections require a reference potential "
          << "if multiple steps.");
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

  // Error check TrialMorph
  if (trial->class_name() == "TrialMorph" || trial->class_name() == "TrialRXNAVBHalf") {
    const int conf = trial->stage(0).select().configuration_index();
    DEBUG("conf  " << conf);
    for (int stage = 0; stage < trial->num_stages(); ++stage) {
      const TrialStage& ts = trial->stage(stage);
      if (ts.select().is_particle_type_set()) {
        const int ptype = ts.select().particle_type();
        DEBUG("ptype " << ptype);
        const int num_sites = system_->configuration(conf).particle_type(ptype).num_sites();
        DEBUG("num_sites " << num_sites);
        if (num_sites > 1) {
          DEBUG("ref " << ts.reference());
          ASSERT(ts.reference() != -1 || !ts.ref().empty(),
            "TrialMorph and TrialRXNAVB requires a reference_index argument when "
            << "the number of sites:" << num_sites << " > 1");
        }
      }
    }
  }

//  // Need to implement some way to handle profiles when config is only in the first select
//  if (trial->num_stages() > 1) {
//    ASSERT(system_->num_configurations() == 1,
//      "not implemented. Fix TrialStage::set_rosenbluth_energy_");
//  }

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
    trial->precompute(criteria_.get(), system_.get());
    trial_factory_->add(trial);
  }

  // HWH deprecate?
  // If later, perhaps after some initialization, more trials are added,
  // then Analyze and Modify classes may need to be re-initialized.
  // analyze_factory_->initialize(criteria_, system_, trial_factory_);
  // modify_factory_->initialize(criteria_, &system_, trial_factory_.get());
}

void MonteCarlo::add(std::shared_ptr<TrialFactoryNamed> trials) {
  ASSERT(criteria_set_, "set Criteria before Trials.");
//  double total_weight = 0.;
//  for (std::shared_ptr<Trial> itrl : trials->trials()) {
//    total_weight += itrl->weight();
//  }
//  DEBUG("total weight " << total_weight);
  for (std::shared_ptr<Trial> itrl : trials->trials()) {
//    DEBUG("itrl weight " << itrl->weight());
//    itrl->set_weight(itrl->weight()*itrl->weight()/total_weight);
//    DEBUG(itrl->weight() << " " << total_weight);
    add(itrl);
  }
}

bool MonteCarlo::duplicate_stepper_output_file_(const std::string output_file) {
  if (!output_file.empty()) {
    for (const std::shared_ptr<Analyze>& an : analyze_factory_->analyzers()) {
      if (an->output_file() == output_file) return true;
    }
    for (const std::shared_ptr<Modify>& mod : modify_factory_->modifiers()) {
      if (mod->output_file() == output_file) return true;
    }
  }
  return false;
}

void MonteCarlo::add(std::shared_ptr<Analyze> analyze) {
  ASSERT(criteria_set_, "set Criteria before Analyze");
  ASSERT(!duplicate_stepper_output_file_(analyze->output_file()),
    "Analyze " << analyze->class_name() << " should not have the same " <<
    "output_file as an already existing Stepper file name: " <<
    analyze->output_file());

  // process multistate
  DEBUG("class name? " << analyze->class_name());
  if (analyze->is_multistate() && analyze->class_name() != "AnalyzeFactory") {
    int trials_per_write = 1;
    std::string output_file;
    if (analyze->is_multistate_aggregate()) {
      trials_per_write = analyze->trials_per_write();
      output_file = analyze->output_file();
      analyze->initialize(this);
    }
    auto multi = MakeAnalyzeFactory({
      {"multistate", "true"},
      {"trials_per_write", str(trials_per_write)},
      {"output_file", output_file},
      {"append", str(analyze->append())},
      {"trials_per_update", "1"},
      {"stop_after_phase", str(analyze->stop_after_phase())},
      {"start_after_phase", str(analyze->start_after_phase())},
      {"stop_after_cycle", str(analyze->stop_after_cycle())},
      {"start_after_cycle", str(analyze->start_after_cycle())},
      {"output_file_append_phase", str(analyze->output_file_append_phase())},
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
        an->empty_output_file();
      }
      an->set_state(state);
      // an->initialize(criteria_, system_, trial_factory_);
      multi->add(an);
    }
    analyze = multi;
  }
  analyze->initialize(this);
  DEBUG("mults " << analyze->is_multistate() << " class name? " << analyze->class_name());
  analyze_factory_->add(analyze);
}

// copied from above
void MonteCarlo::add(std::shared_ptr<Modify> modify) {
  ASSERT(criteria_set_, "set Criteria before Modify");
  ASSERT(!duplicate_stepper_output_file_(modify->output_file()),
    "Modify " << modify->class_name() << " should not have the same " <<
    "output_file as an already existing Stepper file name: " <<
    modify->output_file());
  DEBUG("class name? " << modify->class_name());
  if (modify->is_multistate() && modify->class_name() != "ModifyFactory") {
    int trials_per_write = 1;
    std::string output_file;
    DEBUG("aggregate? " << modify->is_multistate_aggregate());
    if (modify->is_multistate_aggregate()) {
      trials_per_write = modify->trials_per_write();
      output_file = modify->output_file();
      modify->initialize(this);
    }
    auto multi = MakeModifyFactory({
      {"multistate", "true"},
      {"trials_per_write", str(trials_per_write)},
      {"output_file", output_file},
      {"trials_per_update", "1"},
      {"stop_after_phase", str(modify->stop_after_phase())},
      {"start_after_phase", str(modify->start_after_phase())},
      {"stop_after_cycle", str(modify->stop_after_cycle())},
      {"start_after_cycle", str(modify->start_after_cycle())},
      {"output_file_append_phase", str(modify->output_file_append_phase())},
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
        mod->empty_output_file();
      }
      mod->set_state(state);
      // mod->initialize(criteria_.get(), &system_, trial_factory_.get());
      multi->add(mod);
    }
    modify = multi;
  }
  modify->initialize(this);
  DEBUG("mults " << modify->is_multistate() << " class name? " << modify->class_name());

  // Check that modifiers aren't added after ReadConfigFromFile.
  if (modify_factory_->num() > 0) {
    ASSERT(
      modify_factory_->modifiers().back()->class_name() != "ReadConfigFromFile",
      "ReadConfigFromFile should be the last modifier.");
  }
  modify_factory_->add(modify);
}

//void MonteCarlo::add(const std::shared_ptr<Modify> modify) {
//  ASSERT(criteria_set_, "set Criteria before Modify");
//  modify->initialize(criteria_.get(), &system_, trial_factory_.get());
//  modify_factory_->add(modify);
//}

void MonteCarlo::set(const std::shared_ptr<Checkpoint> checkpoint) {
  checkpoint_ = checkpoint;
}

void MonteCarlo::after_trial_modify_() {
  modify_factory_->trial(this);
}

void MonteCarlo::after_trial_checkpoint_() {
  if (checkpoint_) {
    checkpoint_->check(*this);
  }
}

void MonteCarlo::serialize(std::ostream& ostr) const {
  feasst_serialize_version(531, ostr);
  feasst_serialize(system_, ostr);
  feasst_serialize_fstdr(criteria_, ostr);
  feasst_serialize(trial_factory_, ostr);
  feasst_serialize(analyze_factory_, ostr);
  feasst_serialize(modify_factory_, ostr);
  feasst_serialize(checkpoint_, ostr);
  feasst_serialize_fstdr(random_, ostr);
  feasst_serialize_fstdr(action_, ostr);
  feasst_serialize(args_, ostr);
  //feasst_serialize(next_arg_, ostr);
  feasst_serialize(config_set_, ostr);
  feasst_serialize(potential_set_, ostr);
  feasst_serialize(thermo_params_set_, ostr);
  feasst_serialize(system_set_, ostr);
  feasst_serialize(criteria_set_, ostr);
  feasst_serialize(parse_for_num_configs_, ostr);
  feasst_serialize(parse_replace_, ostr);
  feasst_serialize(replace_with_index_, ostr);
  feasst_serialize(timer_, ostr);
  feasst_serialize_endcap("MonteCarlo", ostr);
  DEBUG("size: " << ostr.tellp());
}

MonteCarlo::MonteCarlo(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 529 && version <= 532, "version: " << version);
  feasst_deserialize(system_, istr);
  // feasst_deserialize_fstdr(criteria_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      criteria_ = criteria_->deserialize(istr);
    }
  }
  feasst_deserialize(trial_factory_, istr);
  feasst_deserialize(analyze_factory_, istr);
  feasst_deserialize(modify_factory_, istr);
  // HWH: set check energy every trial for testing: modify_factory_->get_modify(0)->set_trials_per_update(1);
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
  //if (version >= 530) {
  //  feasst_deserialize(&next_arg_, istr);
  //}
  feasst_deserialize(&config_set_, istr);
  feasst_deserialize(&potential_set_, istr);
  feasst_deserialize(&thermo_params_set_, istr);
  feasst_deserialize(&system_set_, istr);
  feasst_deserialize(&criteria_set_, istr);
  if (version >= 530) {
    feasst_deserialize(&parse_for_num_configs_, istr);
    feasst_deserialize(&parse_replace_, istr);
  }
  if (version >= 531) {
    feasst_deserialize(&replace_with_index_, istr);
  }
  feasst_deserialize(timer_, istr);
  feasst_deserialize_endcap("MonteCarlo", istr);
}

void MonteCarlo::load_cache_(const bool load) {
  random_->set_cache_to_load(load);
  system_->load_cache(load);
}

void MonteCarlo::unload_cache_(const MonteCarlo& mc) {
  random_->set_cache_to_unload((*mc.random_));
  system_->unload_cache(mc.system());
}

void MonteCarlo::before_attempts_() {
  ASSERT(system_set_, "system must be set before attempting trials.");
  ASSERT(criteria_set_, "criteria must be set before attempting trials.");
}

void MonteCarlo::revert_(const int trial_index,
    const bool accepted,
    const bool endpoint,
    const bool auto_reject,
    const double ln_prob) {
  trial_factory_->revert(trial_index, accepted, auto_reject, system_.get(), criteria_.get());
  DEBUG("reverting " << criteria_->current_energy());
  criteria_->revert_(accepted, endpoint, ln_prob, trial_factory_->trial(trial_index).accept().updtd());
}

void MonteCarlo::attempt_(int num_trials,
    TrialFactory * trial_factory,
    Random * random) {
  if (trial_factory->num() == 0 && trial_factory->num_attempts() == 0) {
    WARN("No Trials to attempt.");
  }
  before_attempts_();
  for (int trial = 0; trial < num_trials; ++trial) {
    if (timer_) timer_->start(0);
    DEBUG("mc trial: " << trial);
    trial_factory->attempt(criteria_.get(), system_.get(), random);
    if (timer_) timer_->start(1);
    after_trial_analyze_();
    if (timer_) timer_->start(2);
    after_trial_modify_();
    if (timer_) timer_->start(3);
    after_trial_checkpoint_();
  }
  if (timer_) timer_->start(4);
}

bool MonteCarlo::attempt_trial(const int index) {
  return trial_factory_->attempt(criteria_.get(), system_.get(),
                                index, random_.get());
}

void MonteCarlo::imitate_trial_rejection_(const int trial_index,
    const double ln_prob,
    const bool endpoint,
    const bool auto_reject,
    const int state_old,
    const int state_new) {
  trial_factory_->imitate_trial_rejection_(trial_index, auto_reject);
  criteria_->imitate_trial_rejection_(ln_prob, state_old, state_new, endpoint);
}

double MonteCarlo::initialize_system(const int config) {
  return system_->initialize(config);
}

void MonteCarlo::initialize_criteria() {
  ASSERT(criteria_, "err");
  criteria_->initialize(system_.get());
}

void MonteCarlo::initialize_trials() {
  for (int trial = 0; trial < trial_factory_->num(); ++trial) {
    trial_factory_->get_trial(trial)->precompute(criteria_.get(), system_.get());
  }
}

void MonteCarlo::initialize_analyzers() {
  for (int an = 0; an < analyze_factory_->num(); ++an) {
    analyze_factory_->get_analyze(an)->initialize(this);
  }
}

void MonteCarlo::write_checkpoint() const {
  if (checkpoint_) checkpoint_->write(*this);
}

void MonteCarlo::run_until_complete_(TrialFactory * trial_factory,
                                     Random * random) {
  while (!criteria_->is_complete()) {
    DEBUG("here");
    attempt_(1, trial_factory, random);
  }
  write_checkpoint();
  write_to_file();
}

void MonteCarlo::run_until_file_exists(const std::string& file_name,
    const int trials_per_file_check) {
  if (!file_name.empty()) {
    while (!file_exists(file_name)) {
      DEBUG("here");
      attempt_(trials_per_file_check, trial_factory_.get(), random_.get());
    }
    write_checkpoint();
    write_to_file();
  }
}

void MonteCarlo::synchronize_(const MonteCarlo& mc,
     const std::vector<std::shared_ptr<Select> >& perturbed) {
  system_->synchronize_(mc.system(), perturbed);
  criteria_->synchronize_(mc.criteria());
  trial_factory_->synchronize_(mc.trials());
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

void MonteCarlo::adjust_bounds(const bool left_most, const bool right_most,
  const bool left_complete, const bool right_complete,
  const bool all_min_size,
  const int min_size, MonteCarlo * mc) {
  bool adjusted_up;
  std::vector<int> states;
  if (mc) {
    criteria_->adjust_bounds(left_most, right_most, left_complete, right_complete,
      all_min_size, min_size,
      *system_, &mc->system(),  mc->get_criteria(), &adjusted_up, &states);
    DEBUG("adjusted_up " << adjusted_up);
    DEBUG("states: " << feasst_str(states));
    analyze_factory_->adjust_bounds(adjusted_up, states, mc->get_analyze_factory());
    modify_factory_->adjust_bounds(adjusted_up, states, mc->get_modify_factory());
  } else {
    // single processor adjustment on the left and right most only.
    criteria_->adjust_bounds(left_most, right_most, false, false, false, min_size,
      *system_, NULL, NULL, NULL, NULL);
  }
}

void MonteCarlo::ghost_trial_(
    const double ln_prob,
    const int state_old,
    const int state_new,
    const bool endpoint) {
  criteria_->imitate_trial_rejection_(ln_prob, state_old, state_new, endpoint);
}

void MonteCarlo::write_to_file() {
  analyze_factory_->write_to_file(*this);
  modify_factory_->write_to_file(this);
}

void MonteCarlo::run_num_trials(int64_t num_trials) {
  while (num_trials > 0) {
    attempt(1);
    --num_trials;
    DEBUG("num_trials " << num_trials);
  }
}

void MonteCarlo::run_until_num_particles(const int until_num_particles,
    const std::string& particle_type_name, const int configuration_index) {
  const Configuration& conf = configuration(configuration_index);
  int particle_type = -1;
  if (!particle_type_name.empty()) {
    particle_type = conf.particle_name_to_type(particle_type_name);
  }
  while ((until_num_particles > 0) &&
         ((particle_type == -1 && (conf.num_particles() != until_num_particles)) ||
          (particle_type != -1 && (conf.num_particles_of_type(particle_type) != until_num_particles)))) {
    attempt(1);
    DEBUG("num_particles " << conf.num_particles());
  }
}

void MonteCarlo::run_for_hours(const double hours) {
  if (hours > 0) {
    const double begin = cpu_hours();
    while (hours > cpu_hours() - begin) {
      attempt(10);
    }
  }
}

Criteria * MonteCarlo::get_criteria() { return criteria_.get(); }

const Criteria& MonteCarlo::criteria() const {
  ASSERT(criteria_, "Criteria not set.");
  return const_cast<Criteria&>(*criteria_);
}

void MonteCarlo::after_trial_analyze_() {
  analyze_factory_->trial(*this);
}

void MonteCarlo::finalize_(const int trial_index) {
  trial_factory_->finalize(trial_index, system_.get(), criteria_.get());
}

std::string MonteCarlo::serialize() const {
  std::stringstream ss;
  serialize(ss);
  return ss.str();
}
//MonteCarlo MonteCarlo::deserialize(const std::string str) {
//  std::stringstream ss(str);
//  return MonteCarlo(ss);
//}

const Configuration& MonteCarlo::configuration(const int index) const {
  return system_->configuration(index); }
void MonteCarlo::add_to_optimized(std::shared_ptr<Potential> potential,
    const int config) {
  system_->add_to_optimized(potential, config); }
void MonteCarlo::add_to_reference(std::shared_ptr<Potential> potential,
    const int index, const std::string name) {
  system_->add_to_reference(potential, index, name); }
void MonteCarlo::add(std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const int config) {
  system_->add(neighbor_criteria, config); }
const ThermoParams& MonteCarlo::thermo_params() const {
  return system_->thermo_params(); }
const System& MonteCarlo::system() const { return *system_; }
System * MonteCarlo::get_system() { return system_.get(); }
TrialFactory * MonteCarlo::get_trial_factory() { return trial_factory_.get(); }
void MonteCarlo::remove_trial(const int index) { trial_factory_->remove(index); }
const TrialFactory& MonteCarlo::trials() const { return *trial_factory_; }
const Trial& MonteCarlo::trial(const int index) const {
  return trial_factory_->trial(index); }
void MonteCarlo::attempt(const int num_trials) {
  attempt_(num_trials, trial_factory_.get(), random_.get()); }
void MonteCarlo::reset_trial_stats() { trial_factory_->reset_stats(); }
void MonteCarlo::run_until_complete() {
  run_until_complete_(trial_factory_.get(), random_.get()); }
void MonteCarlo::delay_finalize_() {
  trial_factory_->delay_finalize(); }
AnalyzeFactory * MonteCarlo::get_analyze_factory() { return analyze_factory_.get(); }
ModifyFactory * MonteCarlo::get_modify_factory() { return modify_factory_.get(); }
void MonteCarlo::remove_modify(const int index) { modify_factory_->remove(index); }
void MonteCarlo::remove_analyze(const int index) { analyze_factory_->remove(index); }
const Modify& MonteCarlo::modify(const int index) const {
  return modify_factory_->modify(index); }
int MonteCarlo::num_modifiers() const {
  return static_cast<int>(modify_factory_->modifiers().size()); }
const std::vector<std::shared_ptr<Analyze> >& MonteCarlo::analyzers() const {
  return analyze_factory_->analyzers(); }
const Analyze& MonteCarlo::analyze(const int index) const {
  return analyze_factory_->analyze(index); }
int MonteCarlo::num_analyzers() const {
  return static_cast<int>(analyze_factory_->analyzers().size()); }
const std::vector<std::shared_ptr<Modify> >& MonteCarlo::modifiers() const {
  return modify_factory_->modifiers(); }

void MonteCarlo::set_cycles_to_complete(const int num) {
  criteria_->set_cycles_to_complete(num);
}

void MonteCarlo::set_parse_for_num_configs(const int num) {
  DEBUG("num " << num);
  parse_for_num_configs_ = num;
  DEBUG("parse_for_num_configs_ " << parse_for_num_configs_);
}
void MonteCarlo::set_parse_replace(
    const std::vector<std::vector<std::string> >& replace) {
  parse_replace_ = replace;
}
void MonteCarlo::set_replace_with_index(const std::string& str) {
  replace_with_index_ = str;
}

}  // namespace feasst
