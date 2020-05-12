#include <cmath>
#include "monte_carlo/include/trial_select.h"
#include "utils/include/serialize.h"
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
  const TrialSelect& select = (*stages)[0]->trial_select();
  if ((*stages)[0]->rosenbluth().num() > 1) {
    system->get_configuration()->revive(select.mobile());
  }
  acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  DEBUG("old en " << criteria->current_energy());
  DEBUG("new en " << MAX_PRECISION << acceptance->energy_new());
  { // Metropolis
    const Configuration& config = system->configuration();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    DEBUG("selprob " << select.probability() << " betamu " << criteria->beta_mu(particle_type));
    DEBUG("lnselprob " << std::log(select.probability()));
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    acceptance->add_to_ln_metropolis_prob(
      std::log(select.probability())
      + criteria->beta_mu(particle_type)
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
  ASSERT(834 == version, "mismatch version: " << version);
}

void ComputeAddAVB::serialize_compute_add_avb_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(834, ostr);
}

void ComputeAddAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_add_avb_(ostr);
}

}  // namespace feasst
