#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "configuration/include/neighbor_criteria.h"
#include "system/include/energy_map.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_select.h"
#include "cluster/include/compute_add_avb_divalent.h"

namespace feasst {

ComputeAddAVBDivalent::ComputeAddAVBDivalent(argtype args) {
  class_name_ = "ComputeAddAVBDivalent";
  neighbor_ = integer("neighbor_index", &args, 0);
  feasst_check_all_used(args);
}

FEASST_MAPPER(ComputeAddAVBDivalent,);

void ComputeAddAVBDivalent::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeAddAVBDivalent");
  for (int config = 0; config < system->num_configurations(); ++config) {
    ASSERT(acceptance->perturbed(config).num_sites() == 3,
      "hard coded for single site particles");
  }

  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  if (!acceptance->reject()) {
    ASSERT(system->num_configurations() == 1, "not implemented for multiple configs");
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    acceptance->add_to_energy_new(criteria->current_energy());
    //acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
    acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
  //  acceptance->add_to_macrostate_shift(1);
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    DEBUG("old en " << criteria->current_energy());
    DEBUG("new en " << MAX_PRECISION << acceptance->energy_new());

    const Configuration& config = system->configuration();
    const TrialSelect& select0 = (*stages)[0]->trial_select();
    const TrialSelect& select1 = (*stages)[1]->trial_select();
    const TrialSelect& select2 = (*stages)[2]->trial_select();
    const int particle_index0 = select0.mobile().particle_index(0);
    const int particle_type0 = config.select_particle(particle_index0).type();
    const int particle_index1 = select1.mobile().particle_index(0);
    const int particle_type1 = config.select_particle(particle_index1).type();
    const int particle_index2 = select2.mobile().particle_index(0);
    const int particle_type2 = config.select_particle(particle_index2).type();
    ASSERT(particle_type1 == particle_type2,
      "acceptance hard-coded for particle_type1: " << particle_type1 <<
      " == particle_type2: " << particle_type2);
    ASSERT(particle_type0 != particle_type1,
      "acceptance hard-coded for particle_type0: " << particle_type0 <<
      " != particle_type1: " << particle_type1);
    // HWH update with configuration_index_
    const NeighborCriteria& neighbor = system->neighbor_criteria(neighbor_, 0);
    const double volume_av = neighbor.volume(config.dimension());

    select0.map_(*system, neighbor_).neighbors(
      neighbor,
      config,
      select0.mobile().particle_index(0),
      select0.mobile().site_index(0, 0),
      select1.mobile().site_index(0, 0),
      &neighbors_,
      1);
    const int num_neigh = static_cast<int>(neighbors_.num_sites());
    ASSERT(num_neigh > 1, "AVBDivalent shouldn't have one or less neighbors." <<
      " Was a neighbor skipped, for example, using VisitModel::energy_cutoff?");
    const double volume = config.domain().volume();
    const ThermoParams& params = system->thermo_params();
    DEBUG("volume " << volume);
    DEBUG("volume_av " << volume_av);
    DEBUG("num_neigh " << num_neigh);
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    double prob = 1./config.num_particles_of_type(particle_type0);
    DEBUG("prob " << prob);
    DEBUG("nt " << config.num_particles_of_type(particle_type0));
    acceptance->add_to_ln_metropolis_prob(
      std::log(volume*prob)
      + std::log(volume_av/static_cast<double>(num_neigh))
      + std::log(volume_av/static_cast<double>(num_neigh - 1))
      + params.beta_mu(particle_type0)
      + 2*params.beta_mu(particle_type1)
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> ComputeAddAVBDivalent::create(std::istream& istr) const {
  return std::make_shared<ComputeAddAVBDivalent>(istr);
}

ComputeAddAVBDivalent::ComputeAddAVBDivalent(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeAddAVBDivalent", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(5977 == version, "mismatch version: " << version);
  feasst_deserialize(&neighbor_, istr);
}

void ComputeAddAVBDivalent::serialize_compute_add_avb_divalent_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(5977, ostr);
  feasst_serialize(neighbor_, ostr);
}

void ComputeAddAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_add_avb_divalent_(ostr);
}

}  // namespace feasst
