#include "utils/include/serialize.h"
#include "beta_expanded/include/compute_beta.h"

namespace feasst {

ComputeBeta::ComputeBeta() { class_name_ = "ComputeBeta"; }

class MapComputeBeta {
 public:
  MapComputeBeta() {
    auto obj = MakeComputeBeta();
    obj->deserialize_map()["ComputeBeta"] = obj;
  }
};

static MapComputeBeta mapper_ = MapComputeBeta();

void ComputeBeta::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeBeta");
  const double beta_old = system->thermo_params().beta();
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  const double beta_new = system->thermo_params().beta();
  acceptance->set_energy_new(criteria->current_energy());
  acceptance->set_energy_profile_new(criteria->current_energy_profile());
  acceptance->add_to_ln_metropolis_prob(
    -(beta_new - beta_old)*criteria->current_energy());
}

std::shared_ptr<TrialCompute> ComputeBeta::create(std::istream& istr) const {
  return std::make_shared<ComputeBeta>(istr);
}

ComputeBeta::ComputeBeta(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeBeta", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2057 == version, "mismatch version: " << version);
}

void ComputeBeta::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(2057, ostr);
}

}  // namespace feasst
