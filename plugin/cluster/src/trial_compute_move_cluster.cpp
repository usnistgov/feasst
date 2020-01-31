#include "cluster/include/trial_compute_move_cluster.h"

namespace feasst {

TrialComputeMoveCluster::TrialComputeMoveCluster() {
  class_name_ = "TrialComputeMoveCluster";
}

class MapTrialComputeMoveCluster {
 public:
  MapTrialComputeMoveCluster() {
    auto obj = MakeTrialComputeMoveCluster();
    obj->deserialize_map()["TrialComputeMoveCluster"] = obj;
  }
};

static MapTrialComputeMoveCluster mapper_ = MapTrialComputeMoveCluster();

void TrialComputeMoveCluster::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeMoveCluster");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  ASSERT(stages->size() == 1, "assumes 1 stage");
  const TrialSelect * csel = (*stages)[0]->trial_select();
  const int size_old = csel->mobile().num_particles();
  DEBUG("size_old " << size_old);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);

  /** Compute new cluster and require the same cluster size to ensure detailed
      balance. Cluster sizes can only increase as they are rigidly translated
      or rotated. */
  if (!cselect_) {
    std::stringstream ss;
    csel->serialize(ss);
    cselect_ = std::make_shared<TrialSelectCluster>(ss);
  }
  cselect_->select_cluster(csel->mobile().particle_index(0), *system);
  const int size_new = cselect_->mobile().num_particles();
  DEBUG("size_new " << size_new);
  if (size_old != size_new) {
    ASSERT(size_old < size_new, "old size: " << size_old << " must be smaller "
      << "than new size: " << size_new << " unless rigid translate/rotation "
      << "improperly broke apart the cluster");
    acceptance->set_reject(true);
  }
  //ASSERT(size_old == size_new, "detailed balance");

  DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
  DEBUG("new: " << acceptance->energy_new());
  DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
  const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);
}

std::shared_ptr<TrialCompute> TrialComputeMoveCluster::create(std::istream& istr) const {
  return std::make_shared<TrialComputeMoveCluster>(istr);
}

TrialComputeMoveCluster::TrialComputeMoveCluster(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(888 == version, "mismatch version: " << version);
}

void TrialComputeMoveCluster::serialize_trial_compute_move_cluster_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(888, ostr);
}

void TrialComputeMoveCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_move_cluster_(ostr);
}

}  // namespace feasst
