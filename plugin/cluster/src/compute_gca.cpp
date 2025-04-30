#include "utils/include/serialize.h"
#include "monte_carlo/include/acceptance.h"
#include "cluster/include/compute_gca.h"

namespace feasst {

ComputeGCA::ComputeGCA() {
  class_name_ = "ComputeGCA";
}

FEASST_MAPPER(ComputeGCA,);

void ComputeGCA::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeGCA");
  ASSERT(static_cast<int>(stages->size()) == 1, "GCA does not support " <<
    "multiple stages");

  // list pivots, rejected particles and energy of interaction of rejected particles
  // update stage selection and recursively 'compute rosenblut'
  // determine energy change from rejection list
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  ASSERT(system->num_configurations() == 1, "not implemented for multiple configs");

  DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
  DEBUG("new: " << acceptance->energy_new());
  DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());

  const double delta_energy = 0.;
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile());

  // always accept (maybe add some check or force acceptance?)
}

std::shared_ptr<TrialCompute> ComputeGCA::create(std::istream& istr) const {
  return std::make_shared<ComputeGCA>(istr);
}

ComputeGCA::ComputeGCA(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(823 == version, "mismatch version: " << version);
}

void ComputeGCA::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(823, ostr);
}

}  // namespace feasst
