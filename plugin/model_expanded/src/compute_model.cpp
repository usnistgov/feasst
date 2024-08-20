#include "utils/include/serialize.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/acceptance.h"
#include "model_expanded/include/compute_model.h"

namespace feasst {

ComputeModel::ComputeModel() { class_name_ = "ComputeModel"; }

class MapComputeModel {
 public:
  MapComputeModel() {
    auto obj = MakeComputeModel();
    obj->deserialize_map()["ComputeModel"] = obj;
  }
};

static MapComputeModel mapper_ = MapComputeModel();

void ComputeModel::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeModel");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  acceptance->set_energy_new(acceptance->energy_new());
  acceptance->set_energy_profile_new(acceptance->energy_profile_new());
  DEBUG("new energy " << acceptance->energy_new());
  const ThermoParams& thermo = system->thermo_params();
  acceptance->add_to_ln_metropolis_prob(
    // manually add the energy of the old configuration
    // this is an optimization to avoid recomputing the old energy
    + thermo.beta()*criteria->current_energy()
  );
}

std::shared_ptr<TrialCompute> ComputeModel::create(std::istream& istr) const {
  return std::make_shared<ComputeModel>(istr);
}

ComputeModel::ComputeModel(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeModel", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3524 == version, "mismatch version: " << version);
}

void ComputeModel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(3524, ostr);
}

}  // namespace feasst
