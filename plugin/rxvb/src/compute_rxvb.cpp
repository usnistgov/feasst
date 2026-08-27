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

ComputeRXVB::ComputeRXVB(argtype args) : ComputeRXVB(&args) {
  feasst_check_all_used(args);
}
ComputeRXVB::ComputeRXVB(argtype * args) : TrialComputeMove(args) {
  class_name_ = "ComputeRXVB";
  p_bias_ = dble("probability_avb", args, 0.5);
  ASSERT((p_bias_ > 0.) && (p_bias_ < 1.), "probability_avb: " << p_bias_
    << "must be >0 and <1");
  target_particle_type_morph_name_ = str("target_particle_type_morph", args, "-2");
  particle_type_morph_name_ = str("particle_type_morph", args, "-2");
  target_particle_type_name_ = str("target_particle_type", args, "-2");
  particle_type_name_ = str("particle_type", args, "-2");
  avb_ = boolean("avb", args, true);
}

FEASST_MAPPER(ComputeRXVB,);

int ComputeRXVB::num_particles_of_type_(const TrialSelect& select, const Configuration& config) {
  const int particle_index = select.mobile().particle_index(0);
  const int particle_type = config.select_particle(particle_index).type();
  return config.num_particles_of_type(particle_type);
}

void ComputeRXVB::perturb_and_acceptance(
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
    double factor = p_bias_/(1 - p_bias_);
    const int num_targ_part = config.num_particles_of_type(target_particle_type_);
    const int num_targ_part_morph = config.num_particles_of_type(target_particle_type_morph_);
    const int num_part_morph = config.num_particles_of_type(particle_type_morph_);
    DEBUG("num_targ_part " << num_targ_part);
    DEBUG("num_targ_part_morph "<< num_targ_part_morph);
    DEBUG("num_part_morph " << num_part_morph);
    if (avb_) {
      factor = (1 - p_bias_)/p_bias_;
      //ASSERT(std::abs(factor - 1) < NEAR_ZERO, "assumes 1/2 bias");
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
        std::log(factor*num_neighbors*(num_targ_part+1)/num_neigh_new/num_targ_part_morph)
        - system->thermo_params().beta_mu(target_particle_type_)
        - system->thermo_params().beta_mu(particle_type_)
        + system->thermo_params().beta_mu(target_particle_type_morph_)
        + system->thermo_params().beta_mu(particle_type_morph_)
      );
    } else {
      factor = p_bias_/(1 - p_bias_);
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
      //ASSERT(std::abs(factor - 1) < NEAR_ZERO, "assumes 1/2 bias");
      //ASSERT(num_part_morph == 0, "assumes one reactive HWH tmp");
      //ASSERT(num_targ_part_morph == 0, "assumes one reactive HWH tmp");
      acceptance->add_to_ln_metropolis_prob(
        std::log(factor*(num_targ_part_morph+1)*num_neighbors/num_neigh_new/num_targ_part)
        + system->thermo_params().beta_mu(target_particle_type_)
        + system->thermo_params().beta_mu(particle_type_)
        - system->thermo_params().beta_mu(target_particle_type_morph_)
        - system->thermo_params().beta_mu(particle_type_morph_)
      );
    }
  }
}

std::shared_ptr<TrialCompute> ComputeRXVB::create(std::istream& istr) const {
  return std::make_shared<ComputeRXVB>(istr);
}

ComputeRXVB::ComputeRXVB(std::istream& istr)
  : TrialComputeMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(6078 == version, "mismatch version: " << version);
  feasst_deserialize(&p_bias_, istr);
  feasst_deserialize(&avb_, istr);
  feasst_deserialize(&target_particle_type_, istr);
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&target_particle_type_morph_, istr);
  feasst_deserialize(&particle_type_morph_, istr);
  feasst_deserialize(&target_particle_type_name_, istr);
  feasst_deserialize(&particle_type_name_, istr);
  feasst_deserialize(&target_particle_type_morph_name_, istr);
  feasst_deserialize(&particle_type_morph_name_, istr);
}

void ComputeRXVB::serialize_compute_rxnavb_(std::ostream& ostr) const {
  serialize_trial_compute_move_(ostr);
  feasst_serialize_version(6078, ostr);
  feasst_serialize(p_bias_, ostr);
  feasst_serialize(avb_, ostr);
  feasst_serialize(target_particle_type_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(target_particle_type_morph_, ostr);
  feasst_serialize(particle_type_morph_, ostr);
  feasst_serialize(target_particle_type_name_, ostr);
  feasst_serialize(particle_type_name_, ostr);
  feasst_serialize(target_particle_type_morph_name_, ostr);
  feasst_serialize(particle_type_morph_name_, ostr);
}

void ComputeRXVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_rxnavb_(ostr);
}

}  // namespace feasst
