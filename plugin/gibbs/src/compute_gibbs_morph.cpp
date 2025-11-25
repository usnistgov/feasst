#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/perturb.h"
#include "gibbs/include/compute_gibbs_morph.h"

namespace feasst {

ComputeGibbsMorph::ComputeGibbsMorph() {
  class_name_ = "ComputeGibbsMorph";
}

FEASST_MAPPER(ComputeGibbsMorph,);

void ComputeGibbsMorph::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeGibbsMorph");
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());

  // find first and second morph stages
  const TrialStage * first_stage = (*stages)[0];
  const int config_first = first_stage->trial_select().configuration_index();
//  const TrialSelect * select_first = const_cast<const TrialSelect *>(&first_stage->select());
  const int particle_type_first = first_stage->trial_select().particle_type();
//  const TrialSelect * select_second;
  int config_second = -1;
  int particle_type_second = -1;
  bool first = false; // on first PerturbParticleType, reverts to true
  for (TrialStage * stage : *stages) {
    DEBUG("stage name: " << stage->perturb().class_name());
    if (stage->perturb().class_name() == "PerturbParticleType") {
      if (first) {
        config_second = stage->trial_select().configuration_index();
//        select_second = const_cast<const TrialSelect *>(&stage->select());
        particle_type_second = stage->trial_select().particle_type();
      }
      first = !first;
    }
  }
  DEBUG("config_first " << config_first);
  DEBUG("config_second " << config_second);

  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  DEBUG("old energy contribution of first config: " << acceptance->energy_old(config_first));
  DEBUG("old energy contribution of second config: " << acceptance->energy_old(config_second));
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);

  if (!acceptance->reject()) {

    // first
    double delta_energy = acceptance->energy_new(config_first) - acceptance->energy_old(config_first);
    acceptance->set_energy_new(criteria->current_energy(config_first) + delta_energy, config_first);
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config_first), config_first);
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config_first), config_first);
    DEBUG("energy contribution of config " << config_first << ": " << acceptance->energy_old(config_first));

    // second
    delta_energy = acceptance->energy_new(config_second) - acceptance->energy_old(config_second);
    acceptance->set_energy_new(criteria->current_energy(config_second) + delta_energy, config_second);
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config_second), config_second);
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config_second), config_second);
    DEBUG("energy contribution of config " << config_second << ": " << acceptance->energy_new(config_second));

    DEBUG("en 0 current " << MAX_PRECISION << criteria->current_energy(0));
    DEBUG("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
    DEBUG("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));
    DEBUG("en 1 current " << MAX_PRECISION << criteria->current_energy(1));
    DEBUG("en 1 new " << MAX_PRECISION << acceptance->energy_new(1));
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    { // Metropolis
      const Configuration& conf_first = system->configuration(config_first);
      const Configuration& conf_second = system->configuration(config_second);
//      const int particle_first = select_first->mobile().particle_index(0);
//      const int particle_second = select_second->mobile().particle_index(0);
//      DEBUG("particle_first: " << particle_first);
//      DEBUG("particle_second: " << particle_second);
//      // after perturb, types are switched
//      const int particle_type_second = conf_first.select_particle(particle_first).type();
//      const int particle_type_first = conf_second.select_particle(particle_second).type();
//      DEBUG("particle_type_first: " << particle_type_first);
//      DEBUG("particle_type_second: " << particle_type_second);
      const int num_particles_first_in_first = conf_first.num_particles_of_type(particle_type_first);
      const int num_particles_second_in_second = conf_second.num_particles_of_type(particle_type_second);
      const int num_particles_first_in_second = conf_second.num_particles_of_type(particle_type_first);
      const int num_particles_second_in_first = conf_first.num_particles_of_type(particle_type_second);
      DEBUG("num_particles_first_in_first " << num_particles_first_in_first);
      DEBUG("num_particles_second_in_second " << num_particles_second_in_second);
      DEBUG("num_particles_first_in_second " << num_particles_first_in_second);
      DEBUG("num_particles_second_in_first " << num_particles_second_in_first);
      // numbers defined after the morph has already happened
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
      acceptance->add_to_ln_metropolis_prob(
        std::log(static_cast<double>(num_particles_first_in_first + 1)*(num_particles_second_in_second + 1)/num_particles_first_in_second/num_particles_second_in_first)
      );
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    }
  }
}

std::shared_ptr<TrialCompute> ComputeGibbsMorph::create(std::istream& istr) const {
  return std::make_shared<ComputeGibbsMorph>(istr);
}

ComputeGibbsMorph::ComputeGibbsMorph(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeGibbsMorph", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7604 == version, "mismatch version: " << version);
}

void ComputeGibbsMorph::serialize_compute_gibbs_morph_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(7604, ostr);
}

void ComputeGibbsMorph::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_gibbs_morph_(ostr);
}

}  // namespace feasst
