#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_select.h"
#include "cluster/include/compute_move_cluster.h"

namespace feasst {

ComputeMoveCluster::ComputeMoveCluster() {
  class_name_ = "ComputeMoveCluster";
}

class MapComputeMoveCluster {
 public:
  MapComputeMoveCluster() {
    auto obj = MakeComputeMoveCluster();
    obj->deserialize_map()["ComputeMoveCluster"] = obj;
  }
};

static MapComputeMoveCluster mapper_ = MapComputeMoveCluster();

void ComputeMoveCluster::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeMoveCluster");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  ASSERT(stages->size() == 1, "assumes 1 stage");
  const TrialSelect& csel = (*stages)[0]->trial_select();
  const int size_old = csel.mobile().num_particles();
  DEBUG("size_old " << size_old);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);

//  //Use this test to check cluster constraints are working.
//  /** Compute new cluster and require the same cluster size to ensure detailed
//      balance. Cluster sizes can only increase as they are rigidly translated
//      or rotated. */
//  if (!acceptance->reject()) {
//    if (!cselect_) {
//      std::stringstream ss;
//      csel.serialize(ss);
//      cselect_ = std::make_shared<SelectCluster>(ss);
//    }
//    cselect_->select_cluster(csel.mobile().particle_index(0), *system);
//    const int size_new = cselect_->mobile().num_particles();
//    DEBUG("size_new " << size_new);
//    if ((size_old != size_new) && !acceptance->reject()) {
//      DEBUG("acceptance->reject() " << acceptance->reject());
//      ASSERT(size_old < size_new, "old size: " << size_old << " must be smaller "
//        << "than new size: " << size_new << " unless rigid translate/rotation "
//        << "improperly broke apart the cluster");
//      ASSERT(!csel.are_constraints_satisfied(0, *system), "cluster should have changed");
//      DEBUG("reject due to cluster changing");
//      ASSERT(acceptance->reject(), "rej");
//      acceptance->set_reject(true);
//    } else {
//      ASSERT(csel.are_constraints_satisfied(0, *system), "cluster should not have changed");
//      ASSERT(!acceptance->reject(), "rej");
//    }
//    //ASSERT(size_old == size_new, "detailed balance");
//  }

  DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
  DEBUG("new: " << acceptance->energy_new());
  DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
  const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
  acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old());
}

std::shared_ptr<TrialCompute> ComputeMoveCluster::create(std::istream& istr) const {
  return std::make_shared<ComputeMoveCluster>(istr);
}

ComputeMoveCluster::ComputeMoveCluster(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(888 == version, "mismatch version: " << version);
}

void ComputeMoveCluster::serialize_compute_move_cluster_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(888, ostr);
}

void ComputeMoveCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_move_cluster_(ostr);
}

}  // namespace feasst
