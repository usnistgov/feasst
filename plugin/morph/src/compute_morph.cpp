#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/acceptance.h"
#include "morph/include/compute_morph.h"

namespace feasst {

ComputeMorph::ComputeMorph() {
  class_name_ = "ComputeMorph";
}

FEASST_MAPPER(ComputeMorph,);

void ComputeMorph::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeMorph");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  if (!acceptance->reject()) {
    const int config = stages->front()->select().configuration_index();
    DEBUG("old: " << criteria->current_energy(config) << " " << acceptance->energy_old(config));
    DEBUG("new: " << acceptance->energy_new(config));
    DEBUG("energy change: " << acceptance->energy_new(config) - acceptance->energy_old(config));
    const double delta_energy = acceptance->energy_new(config) - acceptance->energy_old(config);
    acceptance->set_energy_new(criteria->current_energy(config) + delta_energy, config);
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config), config);
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config), config);
    const Configuration& configuration = system->configuration(config);

    // initialize delta_
    const int num_ptypes = configuration.num_particle_types();
    if (static_cast<int>(delta_top_.size()) < num_ptypes) {
      delta_top_.resize(num_ptypes);
      delta_bottom_.resize(num_ptypes);
    }
    std::fill(delta_bottom_.begin(), delta_bottom_.end(), 0);

    /// HWH rederive acceptance using number of particles already morphed to/from arrays

    // Metropolis
    const int num_stages = static_cast<int>(stages->size());
    for (int istage = 0; istage < num_stages; ++istage) {
    //for (const TrialStage * stage : *stages) {
      const TrialStage * stage = (*stages)[istage];
      const TrialSelect& select = stage->trial_select();
      const int particle_index = select.mobile().particle_index(0);
      const int particle_type = select.particle_type();
      const int particle_type_morph = configuration.select_particle(particle_index).type();

      ++delta_bottom_[particle_type];
      --delta_bottom_[particle_type_morph];

      std::fill(delta_top_.begin(), delta_top_.end(), 0);
      for (int jstage = istage; jstage < num_stages; ++jstage) {
        const TrialSelect& jselect = (*stages)[jstage]->trial_select();
        const int particle_jindex = jselect.mobile().particle_index(0);
        const int particle_jtype = jselect.particle_type();
        const int particle_jtype_morph = configuration.select_particle(particle_jindex).type();
        ++delta_top_[particle_jtype];
        --delta_top_[particle_jtype_morph];
      }

      const double num_type = configuration.num_particles_of_type(particle_type);
      const double num_type_morph = configuration.num_particles_of_type(particle_type_morph);
      ASSERT(particle_type != particle_type_morph, "err");
      DEBUG("betamu " << system->thermo_params().beta_mu(particle_type));
      acceptance->add_to_ln_metropolis_prob(
        std::log((num_type + delta_top_[particle_type])/(
                  num_type_morph + 1 + delta_bottom_[particle_type_morph]))
        - system->thermo_params().beta_mu(particle_type)
        + system->thermo_params().beta_mu(particle_type_morph)
      );
    }
  }
}

std::shared_ptr<TrialCompute> ComputeMorph::create(std::istream& istr) const {
  return std::make_shared<ComputeMorph>(istr);
}

ComputeMorph::ComputeMorph(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeMorph", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6389 == version, "mismatch version: " << version);
}

void ComputeMorph::serialize_trial_compute_add_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(6389, ostr);
}

void ComputeMorph::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_add_(ostr);
}

}  // namespace feasst
