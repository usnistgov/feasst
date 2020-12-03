#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialComputeAdd::TrialComputeAdd(const argtype& args) : TrialCompute(args) {
  class_name_ = "TrialComputeAdd";
}

class MapTrialComputeAdd {
 public:
  MapTrialComputeAdd() {
    auto obj = MakeTrialComputeAdd();
    obj->deserialize_map()["TrialComputeAdd"] = obj;
  }
};

static MapTrialComputeAdd mapper_ = MapTrialComputeAdd();

void TrialComputeAdd::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeAdd");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();
    const TrialSelect& select = (*stages)[0]->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    DEBUG("volume " << volume << " selprob " << select.probability() << " betamu " << system->thermo_params().beta_mu(particle_type));
    acceptance->add_to_ln_metropolis_prob(
      std::log(volume*select.probability())
      + system->thermo_params().beta_mu(particle_type)
    );
  }
}

std::shared_ptr<TrialCompute> TrialComputeAdd::create(std::istream& istr) const {
  return std::make_shared<TrialComputeAdd>(istr);
}

TrialComputeAdd::TrialComputeAdd(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(834 == version, "mismatch version: " << version);
}

void TrialComputeAdd::serialize_trial_compute_add_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(834, ostr);
}

void TrialComputeAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_add_(ostr);
}

}  // namespace feasst
