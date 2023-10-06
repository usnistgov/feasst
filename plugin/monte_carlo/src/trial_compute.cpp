#include <vector>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

TrialCompute::TrialCompute(argtype args) : TrialCompute(&args) {
  FEASST_CHECK_ALL_USED(args); }
TrialCompute::TrialCompute(argtype * args) {}

void TrialCompute::compute_rosenbluth(
    const int old,
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  const bool is_new_only = (*stages)[0]->is_new_only();
  if (old == 1 && is_new_only) {
    for (TrialStage* stage : *stages) {
      stage->set_mobile_physical(true, system);
//      stage->get_trial_select()->set_trial_state(2);
    }
    return;
  }
  double ln_rosenbluth = 0.;
  double energy_change = 0.;
  bool reference_used = false;
  std::vector<int> configs_used;
  for (TrialStage* stage : *stages) {
    DEBUG("*** Attempting stage. old: " << old << " ***");
    stage->attempt(system, acceptance, criteria, old, random);
    if (stage->rosenbluth().chosen_step() == -1) {
      if (!is_new_only) {
        acceptance->set_reject(true);
        DEBUG("auto reject");
        for (TrialStage* stage : *stages) {
          stage->set_mobile_physical(true, system);
        }
        return;
      }
    }
    if (!stage->are_constraints_satisfied(old, *system)) {
      // reject using ln_prob to incorporate tuning
      acceptance->add_to_ln_metropolis_prob(-1e100);
      //acceptance->set_reject(true);
    }
    double energy;
    const int config = stage->trial_select().configuration_index();
    if (!find_in_list(config, configs_used)) {
      configs_used.push_back(config);
    }
    if (old == 1) {
      energy = stage->rosenbluth().energy(0);
      acceptance->add_to_energy_old(energy, config);
      acceptance->add_to_energy_profile_old(stage->rosenbluth().energy_profile(0), config);
      ln_rosenbluth -= stage->rosenbluth().ln_total_rosenbluth();
      DEBUG("adding to old energy " << energy);
    } else {
      energy = stage->rosenbluth().chosen_energy();
      DEBUG("energy new " << acceptance->energy_new());
      DEBUG("energy " << energy);
      acceptance->add_to_energy_new(energy, config);
      acceptance->add_to_energy_profile_new(stage->rosenbluth().chosen_energy_profile(), config);
      DEBUG("energy new updated " << acceptance->energy_new());
      ln_rosenbluth += stage->rosenbluth().ln_total_rosenbluth();
      DEBUG("adding to new energy " << energy);
    }
    energy_change += energy;
    if (stage->reference() >= 0) reference_used = true;
  }

  // update the trial state of the perturbed selection for each config
  // first, find the first stage for each config in acceptance
  // use each of those first stages to set trial_state
  int trial_state = -1;
  if (acceptance->num_configurations() == 1) {
    trial_state = (*stages)[0]->trial_select().mobile().trial_state();
    // set the trial state if old configuration and is a move type (1)
    if (trial_state == 1 && old == 1) trial_state = 0;
    const int config = (*stages)[0]->trial_select().configuration_index();
    acceptance->set_perturbed_state(trial_state, config);
  } else {
    DEBUG("acceptance->num_configurations() " << acceptance->num_configurations());
    for (int iconf = 0; iconf < acceptance->num_configurations(); ++iconf) {
      DEBUG("iconf " << iconf);
      bool found = false;
      for (TrialStage * stage : *stages) {
        if (!found && iconf == stage->select().configuration_index()) {
          trial_state = stage->select().mobile().trial_state();
          DEBUG("trial state " << trial_state);
          acceptance->set_perturbed_state(trial_state, iconf);
          found = true;
        }
      }
    }
  }

  DEBUG("reference used? " << reference_used);
  if (reference_used) {
    for (int conf = 0; conf < acceptance->num_configurations(); ++conf) {
      if (acceptance->updated(conf) == 1 && find_in_list(conf, configs_used)) {
        DEBUG("conf " << conf);
        //ASSERT(acceptance->perturbed(conf).num_sites() > 0, "error");
        if (acceptance->perturbed(conf).num_sites() > 0) {
          DEBUG(acceptance->perturbed(conf).str());
          DEBUG("state " << acceptance->perturbed(conf).trial_state());
          const double en_full = system->perturbed_energy(acceptance->perturbed(conf), conf);
          const std::vector<double>& en_profile_full = system->stored_energy_profile(conf);
          DEBUG("en_full: " << en_full);
          DEBUG("en_profile_full: " << feasst_str(en_profile_full));
          DEBUG("energy ref: " << energy_change);
          acceptance->set_energy_ref(energy_change, conf);
          DEBUG("old " << old);
          DEBUG("trial_state " << trial_state);
          if (old == 1) {
            acceptance->set_energy_old(en_full, conf);
            acceptance->set_energy_profile_old(en_profile_full, conf);
            const double ln_met = -1.*system->thermo_params().beta()*
              (-en_full + energy_change);
            DEBUG("ln_met " << ln_met);
            acceptance->add_to_ln_metropolis_prob(ln_met);
          } else {
            acceptance->set_energy_new(en_full, conf);
            acceptance->set_energy_profile_new(en_profile_full, conf);
            DEBUG("updated energy_new " << acceptance->energy_new(conf));
            const double ln_met = -1.*system->thermo_params().beta()*
              (en_full - energy_change);
            DEBUG("ln_met " << ln_met);
            acceptance->add_to_ln_metropolis_prob(ln_met);
          }
        }
      }
    }
  }
  DEBUG("ln_rosenbluth " << ln_rosenbluth);
  acceptance->add_to_ln_metropolis_prob(ln_rosenbluth);
}

std::map<std::string, std::shared_ptr<TrialCompute> >& TrialCompute::deserialize_map() {
  static std::map<std::string, std::shared_ptr<TrialCompute> >* ans =
     new std::map<std::string, std::shared_ptr<TrialCompute> >();
  return *ans;
}

void TrialCompute::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<TrialCompute> TrialCompute::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<TrialCompute> TrialCompute::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void TrialCompute::serialize_trial_compute_(std::ostream& ostr) const {
  feasst_serialize_version(803, ostr);
}

TrialCompute::TrialCompute(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(803 == version, "mismatch version: " << version);
}

}  // namespace feasst
