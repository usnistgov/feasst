#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "system/include/neighbor_criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "cluster/include/compute_add_avb_divalent.h"

namespace feasst {

ComputeAddAVBDivalent::ComputeAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria) {
  class_name_ = "ComputeAddAVBDivalent";
  neighbor_criteria_ = neighbor_criteria;
}

class MapComputeAddAVBDivalent {
 public:
  MapComputeAddAVBDivalent() {
    auto obj = MakeComputeAddAVBDivalent(MakeNeighborCriteria());
    obj->deserialize_map()["ComputeAddAVBDivalent"] = obj;
  }
};

static MapComputeAddAVBDivalent mapper_ = MapComputeAddAVBDivalent();

void ComputeAddAVBDivalent::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeAddAVBDivalent");

  ASSERT(acceptance->perturbed().num_sites() == 3,
    "hard coded for single site particles");

  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  DEBUG("old en " << criteria->current_energy());
  DEBUG("new en " << MAX_PRECISION << acceptance->energy_new());

  const Configuration& config = system->configuration();
  const TrialSelect& select0 = (*stages)[0]->trial_select();
  const TrialSelect& select1 = (*stages)[1]->trial_select();
  //const TrialSelect& select2 = (*stages)[2]->trial_select();
  const int particle_index0 = select0.mobile().particle_index(0);
  const int particle_type0 = config.select_particle(particle_index0).type();
  const int particle_index1 = select1.mobile().particle_index(0);
  const int particle_type1 = config.select_particle(particle_index1).type();
  const double volume_av = neighbor_criteria_->volume(config.dimension());

  select0.map_(*system, *neighbor_criteria_).neighbors(
    *neighbor_criteria_,
    config,
    select0.mobile().particle_index(0),
    select0.mobile().site_index(0, 0),
    select1.mobile().site_index(0, 0),
    random,
    &neighbors_,
    1);
  const int num_neigh = static_cast<int>(neighbors_.num_sites());
  const double volume = config.domain().volume();
  //set_probability(volume_av/static_cast<double>(num_neighbors + 1 + delta_ab));
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  acceptance->add_to_ln_metropolis_prob(
    std::log(volume/(config.num_particles_of_type(particle_type0)))
    + std::log(volume_av/static_cast<double>(num_neigh))
    + std::log(volume_av/static_cast<double>(num_neigh - 1))
    + criteria->beta_mu(particle_type0)
    + 2*criteria->beta_mu(particle_type1)
  );
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
}

std::shared_ptr<TrialCompute> ComputeAddAVBDivalent::create(std::istream& istr) const {
  return std::make_shared<ComputeAddAVBDivalent>(istr);
}

ComputeAddAVBDivalent::ComputeAddAVBDivalent(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeAddAVBDivalent", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(5977 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work
  // feasst_deserialize(neighbor_criteria_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      neighbor_criteria_ = std::make_shared<NeighborCriteria>(istr);
    }
  }
}

void ComputeAddAVBDivalent::serialize_compute_add_avb_divalent_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(5977, ostr);
  feasst_serialize(neighbor_criteria_, ostr);
}

void ComputeAddAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_add_avb_divalent_(ostr);
}

}  // namespace feasst
