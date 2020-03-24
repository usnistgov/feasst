#include <cmath>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_compute_remove.h"

namespace feasst {

TrialComputeRemove::TrialComputeRemove() {
  class_name_ = "TrialComputeRemove";
}

class MapTrialComputeRemove {
 public:
  MapTrialComputeRemove() {
    auto obj = MakeTrialComputeRemove();
    obj->deserialize_map()["TrialComputeRemove"] = obj;
  }
};

static MapTrialComputeRemove mapper_ = MapTrialComputeRemove();

void TrialComputeRemove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeRemove");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
  acceptance->add_to_macrostate_shift(-1);
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain()->volume();
    const TrialSelect * select = (*stages)[0]->trial_select();
    const int particle_index = select->mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    acceptance->set_macrostate_shift_type(particle_type);
    DEBUG("volume " << volume << " selprob " << select->probability() << " betamu " << criteria->beta_mu(particle_type));
    acceptance->add_to_ln_metropolis_prob(
      - std::log(volume*select->probability())
      - criteria->beta_mu(particle_type)
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> TrialComputeRemove::create(
    std::istream& istr) const {
  return std::make_shared<TrialComputeRemove>(istr);
}

TrialComputeRemove::TrialComputeRemove(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeRemove", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(332 == version, "mismatch version: " << version);
}

void TrialComputeRemove::serialize_trial_compute_remove_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(332, ostr);
}

void TrialComputeRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_remove_(ostr);
}

}  // namespace feasst
