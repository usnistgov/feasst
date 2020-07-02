#include <vector>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

TrialCompute::TrialCompute() {}

void TrialCompute::compute_rosenbluth(
    const int old,
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  double ln_rosenbluth = 0.;
  double energy_change = 0.;
  bool reference_used = false;
  for (TrialStage* stage : *stages) {
    DEBUG("*** Attempting stage ***");
    stage->attempt(system, acceptance, criteria, old, random);
    if (stage->rosenbluth().chosen_step() == -1) {
      if (!stage->is_new_only()) {
        acceptance->set_reject(true);
        DEBUG("auto reject");
        for (TrialStage* stage : *stages) {
          stage->set_mobile_physical(true, system);
        }
        return;
      }
    }
    if (!stage->are_constraints_satisfied(*system)) {
      acceptance->set_reject(true);
    }
    double energy;
    if (old == 1) {
      energy = stage->rosenbluth().energy(0);
      acceptance->add_to_energy_old(energy);
      ln_rosenbluth -= stage->rosenbluth().ln_total_rosenbluth();
      DEBUG("adding to old energy " << energy);
    } else {
      energy = stage->rosenbluth().chosen_energy();
      acceptance->add_to_energy_new(energy);
      ln_rosenbluth += stage->rosenbluth().ln_total_rosenbluth();
      DEBUG("adding to new energy " << energy);
    }
    energy_change += energy;
    if (stage->reference() >= 0) {
      reference_used = true;
    }
  }
  DEBUG("reference used? " << reference_used);
  if (reference_used) {
    ASSERT(acceptance->perturbed().num_sites() > 0, "error");
    int trial_state = (*stages)[0]->trial_select().mobile().trial_state();
    // set the trial state if old configuration and is a move type (1)
    if (trial_state == 1 && old == 1) trial_state = 0;
    acceptance->set_perturbed_state(trial_state);
    DEBUG(acceptance->perturbed().str());
//    { // delete me
//      for (int p = 0; p < acceptance->perturbed().num_particles(); ++p) {
//        const int pindex = acceptance->perturbed().particle_index(p);
//        for (int s : acceptance->perturbed().site_indices(p)) {
//          DEBUG("ps " << p << " " << s << " ph " << system->configuration().select_particle(pindex).site(s).is_physical());
//        }
//      }
//    }
    DEBUG("state " << acceptance->perturbed().trial_state());
    const double en_full = system->perturbed_energy(acceptance->perturbed());
    DEBUG("en_full: " << en_full);
    DEBUG("energy ref: " << energy_change);
    acceptance->set_energy_ref(energy_change);
    if (old == 1) {
      acceptance->set_energy_old(en_full);
      acceptance->add_to_ln_metropolis_prob(-1.*criteria->beta()*
        (-en_full + energy_change));
    } else {
      acceptance->set_energy_new(en_full);
      acceptance->add_to_ln_metropolis_prob(-1.*criteria->beta()*
        (en_full - energy_change));
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
