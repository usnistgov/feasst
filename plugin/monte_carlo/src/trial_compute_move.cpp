#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

TrialComputeMove::TrialComputeMove() {
  class_name_ = "TrialComputeMove";
}

class MapTrialComputeMove {
 public:
  MapTrialComputeMove() {
    auto obj = MakeTrialComputeMove();
    obj->deserialize_map()["TrialComputeMove"] = obj;
  }
};

static MapTrialComputeMove mapper_ = MapTrialComputeMove();

void TrialComputeMove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeMove");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  for (TrialStage * stage : *stages) stage->mid_stage(system);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
  DEBUG("new: " << acceptance->energy_new());
  DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
  const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);
}

std::shared_ptr<TrialCompute> TrialComputeMove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeMove>(istr);
}

TrialComputeMove::TrialComputeMove(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeMove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(888 == version, "mismatch version: " << version);
}

void TrialComputeMove::serialize_trial_compute_move_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(888, ostr);
}

void TrialComputeMove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_move_(ostr);
}

}  // namespace feasst
