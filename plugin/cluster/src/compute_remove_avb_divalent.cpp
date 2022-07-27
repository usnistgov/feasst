#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "configuration/include/neighbor_criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "cluster/include/compute_remove_avb_divalent.h"

namespace feasst {

ComputeRemoveAVBDivalent::ComputeRemoveAVBDivalent(argtype args) {
  class_name_ = "ComputeRemoveAVBDivalent";
  neighbor_ = integer("neighbor_index", &args, 0);
  FEASST_CHECK_ALL_USED(args);
}

class MapComputeRemoveAVBDivalent {
 public:
  MapComputeRemoveAVBDivalent() {
    auto obj = MakeComputeRemoveAVBDivalent();
    obj->deserialize_map()["ComputeRemoveAVBDivalent"] = obj;
  }
};

static MapComputeRemoveAVBDivalent mapper_ = MapComputeRemoveAVBDivalent();

void ComputeRemoveAVBDivalent::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeRemoveAVBDivalent");

  ASSERT(acceptance->perturbed().num_sites() == 3,
    "hard coded for single site particles");

  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
  acceptance->set_energy_profile_new(criteria->current_energy_profile());
  acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old());
  acceptance->add_to_macrostate_shift(-1);
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

  DEBUG("sel0 " << select0.mobile().str());
  DEBUG("sel1 " << select1.mobile().str());
  DEBUG("sel2 " << select2.mobile().str());

  // HWH Optimize: use select[1,2].probability to obtain number of neighbors.
  Select neighbors_;
  const NeighborCriteria& neighbor = system->neighbor_criteria(neighbor_);
  select0.map_(*system, neighbor_).neighbors(
    neighbor,
    config,
    select0.mobile().particle_index(0),
    select0.mobile().site_index(0, 0),
    select1.mobile().site_index(0, 0),
    &neighbors_);
  const int num_neigh = static_cast<int>(neighbors_.num_sites());
  const double volume_av = neighbor.volume(config.dimension());
//  DEBUG("num_neigh " << num_neigh);
  const double volume = config.domain().volume();
  const ThermoParams& params = system->thermo_params();
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  acceptance->add_to_ln_metropolis_prob(
    - std::log(volume*select0.probability())
    - std::log(volume_av/static_cast<double>(num_neigh))
    - std::log(volume_av/static_cast<double>(num_neigh - 1))
    //-std::log(select1.probability())
    //-std::log(select2.probability())
    - params.beta_mu(particle_type0)
    - 2*params.beta_mu(particle_type1)
  );
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
}

std::shared_ptr<TrialCompute> ComputeRemoveAVBDivalent::create(std::istream& istr) const {
  return std::make_shared<ComputeRemoveAVBDivalent>(istr);
}

ComputeRemoveAVBDivalent::ComputeRemoveAVBDivalent(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeRemoveAVBDivalent", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(5977 == version, "mismatch version: " << version);
  feasst_deserialize(&neighbor_, istr);
}

void ComputeRemoveAVBDivalent::serialize_compute_remove_avb_divalent_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(5977, ostr);
  feasst_serialize(neighbor_, ostr);
}

void ComputeRemoveAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_remove_avb_divalent_(ostr);
}

}  // namespace feasst
