#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
//#include "math/include/constants.h"
#include "configuration/include/select.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "system/include/energy_map.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_select.h"
#include "rxvb/include/compute_rxvb.h"

namespace feasst {

ComputeRxVB::ComputeRxVB(argtype args) : ComputeRxVB(&args) {
  feasst_check_all_used(args);
}
ComputeRxVB::ComputeRxVB(argtype * args) : TrialComputeMove(args) {
  class_name_ = "ComputeRxVB";
  target_particle_type_morph_name_ = str("target_particle_type_morph", args, "-2");
  particle_type_morph_name_ = str("particle_type_morph", args, "-2");
  target_particle_type_name_ = str("target_particle_type", args, "-2");
  particle_type_name_ = str("particle_type", args, "-2");
  if (used("probability_avb", *args)) {
    WARN("ComputeRxVB::probability_avb is deprecated.");
  }
  if (used("avb", *args)) {
    WARN("ComputeRxVB::avb is deprecated. Use forward instead.");
    forward_ = boolean("avb", args);
    ASSERT(!used("forward", *args), "Cannot use both avb and forward.");
  } else {
    forward_ = boolean("forward", args, true);
  }
}

FEASST_MAPPER(ComputeRxVB,);

int ComputeRxVB::num_particles_of_type_(const TrialSelect& select, const Configuration& config) {
  const int particle_index = select.mobile().particle_index(0);
  const int particle_type = config.select_particle(particle_index).type();
  return config.num_particles_of_type(particle_type);
}

void ComputeRxVB::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);

  if (!acceptance->reject()) {
    ASSERT(system->num_configurations() == 1,
      "not implemented for multiple configs");
    const int iconf = stages->front()->select().configuration_index();
    DEBUG("iconf " << iconf);
    const Configuration& config = system->configuration(iconf);

    // initialize particle types
    if (target_particle_type_morph_ == -2) {
      target_particle_type_morph_ = config.particle_name_to_type(target_particle_type_morph_name_);
      particle_type_morph_ = config.particle_name_to_type(particle_type_morph_name_);
      target_particle_type_ = config.particle_name_to_type(target_particle_type_name_);
      particle_type_ = config.particle_name_to_type(particle_type_name_);
    }

    DEBUG("mu("<<target_particle_type_<<")="<<system->thermo_params().beta_mu(target_particle_type_));
    DEBUG("mu("<<particle_type_<<")="<<system->thermo_params().beta_mu(particle_type_));
    DEBUG("mu("<<target_particle_type_morph_<<")="<<system->thermo_params().beta_mu(target_particle_type_morph_));
    DEBUG("mu("<<particle_type_morph_<<")="<<system->thermo_params().beta_mu(particle_type_morph_));

    // obtain delta energy
    DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
    DEBUG("new: " << acceptance->energy_new());
    DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
    const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
    acceptance->set_energy_new(criteria->current_energy() + delta_energy);
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
    acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old());

    // obtain Metropolis prefactors
    const TrialSelect& select0 = (*stages)[0]->trial_select();
    const TrialSelect& select1 = (*stages)[1]->trial_select();
    const int num_targ_part = config.num_particles_of_type(target_particle_type_);
    const int num_targ_part_morph = config.num_particles_of_type(target_particle_type_morph_);
    const int num_part_morph = config.num_particles_of_type(particle_type_morph_);
    DEBUG("num_targ_part " << num_targ_part);
    DEBUG("num_targ_part_morph "<< num_targ_part_morph);
    DEBUG("num_part_morph " << num_part_morph);
    if (forward_) {
      const double num_neighbors = select0.probability();
      DEBUG("num_neighbors " << num_neighbors);
      const int neighbor_index = 0; // assumed
      const NeighborCriteria& neighbor = system->neighbor_criteria(neighbor_index, 0);
      select0.map_(*system, neighbor_index).neighbors(
        neighbor,
        config,
        select1.mobile().particle_index(0),
        select1.mobile().site_index(0, 0),
        select0.mobile().site_index(0, 0),
        &neighbors_,
        1);
      const int num_neigh_new = static_cast<int>(neighbors_.num_sites());
      DEBUG("num_neigh_new " << num_neigh_new);
      ASSERT(num_neigh_new > 0, "There should be atleast one neighbor because we put one there!");
      //ASSERT(num_part_morph == 1, "assumes one reactive HWH tmp");
      //ASSERT(num_targ_part_morph == 1, "assumes one reactive HWH tmp");
      acceptance->add_to_ln_metropolis_prob(
        std::log(num_neighbors*(num_targ_part+1)/num_neigh_new/num_targ_part_morph)
        - system->thermo_params().beta_mu(target_particle_type_)
        - system->thermo_params().beta_mu(particle_type_)
        + system->thermo_params().beta_mu(target_particle_type_morph_)
        + system->thermo_params().beta_mu(particle_type_morph_)
      );
    } else {
      const double num_neighbors = select0.probability();
      DEBUG("num_neighbors " << num_neighbors);
      const int neighbor_index = 0; // assumed
      const NeighborCriteria& neighbor = system->neighbor_criteria(neighbor_index, 0);
      select0.map_(*system, neighbor_index).neighbors(
        neighbor,
        config,
        select1.mobile().particle_index(0),
        select1.mobile().site_index(0, 0),
        select0.mobile().site_index(0, 0),
        &neighbors_,
        1);
      const int num_neigh_new = static_cast<int>(neighbors_.num_sites());
      DEBUG("num_neigh_new " << num_neigh_new);
      ASSERT(num_neigh_new > 0, "There should be atleast one neighbor because we put one there!");
      //ASSERT(num_part_morph == 0, "assumes one reactive HWH tmp");
      //ASSERT(num_targ_part_morph == 0, "assumes one reactive HWH tmp");
      acceptance->add_to_ln_metropolis_prob(
        std::log((num_targ_part_morph+1)*num_neighbors/num_neigh_new/num_targ_part)
        + system->thermo_params().beta_mu(target_particle_type_)
        + system->thermo_params().beta_mu(particle_type_)
        - system->thermo_params().beta_mu(target_particle_type_morph_)
        - system->thermo_params().beta_mu(particle_type_morph_)
      );
    }
  }
}

std::shared_ptr<TrialCompute> ComputeRxVB::create(std::istream& istr) const {
  return std::make_shared<ComputeRxVB>(istr);
}

ComputeRxVB::ComputeRxVB(std::istream& istr)
  : TrialComputeMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 6078 && version <= 6079, "mismatch version: " << version);
  if (version == 6078) {
    double p_bias;
    feasst_deserialize(&p_bias, istr);
  }
  feasst_deserialize(&forward_, istr);
  feasst_deserialize(&target_particle_type_, istr);
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&target_particle_type_morph_, istr);
  feasst_deserialize(&particle_type_morph_, istr);
  feasst_deserialize(&target_particle_type_name_, istr);
  feasst_deserialize(&particle_type_name_, istr);
  feasst_deserialize(&target_particle_type_morph_name_, istr);
  feasst_deserialize(&particle_type_morph_name_, istr);
}

void ComputeRxVB::serialize_compute_rxvb_(std::ostream& ostr) const {
  serialize_trial_compute_move_(ostr);
  feasst_serialize_version(6079, ostr);
  feasst_serialize(forward_, ostr);
  feasst_serialize(target_particle_type_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(target_particle_type_morph_, ostr);
  feasst_serialize(particle_type_morph_, ostr);
  feasst_serialize(target_particle_type_name_, ostr);
  feasst_serialize(particle_type_name_, ostr);
  feasst_serialize(target_particle_type_morph_name_, ostr);
  feasst_serialize(particle_type_morph_name_, ostr);
}

void ComputeRxVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_rxvb_(ostr);
}

}  // namespace feasst
