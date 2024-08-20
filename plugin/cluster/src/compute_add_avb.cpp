#include <cmath>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/acceptance.h"
#include "cluster/include/compute_add_avb.h"

namespace feasst {

ComputeAddAVB::ComputeAddAVB() {
  class_name_ = "ComputeAddAVB";
}

class MapComputeAddAVB {
 public:
  MapComputeAddAVB() {
    auto obj = MakeComputeAddAVB();
    obj->deserialize_map()["ComputeAddAVB"] = obj;
  }
};

static MapComputeAddAVB mapper_ = MapComputeAddAVB();

void ComputeAddAVB::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeAddAVB");
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  acceptance->add_to_energy_new(criteria->current_energy());
  //acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
  acceptance->add_to_macrostate_shift(1);
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  DEBUG("old en " << criteria->current_energy());
  DEBUG("new en " << MAX_PRECISION << acceptance->energy_new());
  { // Metropolis
    const Configuration& config = system->configuration();
    const TrialSelect& select = (*stages)[0]->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    acceptance->set_macrostate_shift_type(particle_type);
    const ThermoParams& params = system->thermo_params();
    DEBUG("selprob " << select.probability() << " betamu " << params.beta_mu(particle_type));
    DEBUG("lnselprob " << std::log(select.probability()));
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    acceptance->add_to_ln_metropolis_prob(
      std::log(select.probability())
      + params.beta_mu(particle_type)
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> ComputeAddAVB::create(std::istream& istr) const {
  return std::make_shared<ComputeAddAVB>(istr);
}

ComputeAddAVB::ComputeAddAVB(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeAddAVB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(1676 == version, "mismatch version: " << version);
}

void ComputeAddAVB::serialize_compute_add_avb_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(1676, ostr);
}

void ComputeAddAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_add_avb_(ostr);
}

}  // namespace feasst
