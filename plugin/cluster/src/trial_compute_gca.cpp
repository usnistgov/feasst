#include "cluster/include/trial_compute_gca.h"
#include "utils/include/serialize.h"

namespace feasst {

TrialComputeGCA::TrialComputeGCA() {
  class_name_ = "TrialComputeGCA";
}

class MapTrialComputeGCA {
 public:
  MapTrialComputeGCA() {
    auto obj = MakeTrialComputeGCA();
    obj->deserialize_map()["TrialComputeGCA"] = obj;
  }
};

static MapTrialComputeGCA mapper_ = MapTrialComputeGCA();

void TrialComputeGCA::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeGCA");
  ASSERT(static_cast<int>(stages->size()) == 1, "GCA does not support " <<
    "multiple stages");

  // list pivots, rejected particles and energy of interaction of rejected particles
  // update stage selection and recursively 'compute rosenblut'
  // determine energy change from rejection list
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);

  DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
  DEBUG("new: " << acceptance->energy_new());
  DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());

  const double delta_energy = 0.;
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);

  // always accept (maybe add some check or force acceptance?)
}

std::shared_ptr<TrialCompute> TrialComputeGCA::create(std::istream& istr) const {
  return std::make_shared<TrialComputeGCA>(istr);
}

TrialComputeGCA::TrialComputeGCA(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(823 == version, "mismatch version: " << version);
}

void TrialComputeGCA::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(823, ostr);
}

}  // namespace feasst
