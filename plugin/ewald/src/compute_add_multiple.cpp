#include <cmath>
#include "monte_carlo/include/trial_select.h"
#include "utils/include/serialize.h"
#include "ewald/include/compute_add_multiple.h"

namespace feasst {

ComputeAddMultiple::ComputeAddMultiple() {
  class_name_ = "ComputeAddMultiple";
}

class MapComputeAddMultiple {
 public:
  MapComputeAddMultiple() {
    auto obj = MakeComputeAddMultiple();
    obj->deserialize_map()["ComputeAddMultiple"] = obj;
  }
};

static MapComputeAddMultiple mapper_ = MapComputeAddMultiple();

void ComputeAddMultiple::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeAddMultiple");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  DEBUG("deltaE " << MAX_PRECISION << acceptance->energy_new());
  for (const TrialStage * stage : *stages) {
    { // Metropolis
      const Configuration& config = system->configuration();
      const double volume = config.domain().volume();
      const TrialSelect& select = stage->trial_select();
      const int particle_index = select.mobile().particle_index(0);
      const int particle_type = config.select_particle(particle_index).type();
      DEBUG("volume " << volume << " selprob " << select.probability() <<
        " betamu " << criteria->beta_mu(particle_type));
      acceptance->add_to_ln_metropolis_prob(
        std::log(volume*select.probability())
        + criteria->beta_mu(particle_type)
      );
    }
  }
}

std::shared_ptr<TrialCompute> ComputeAddMultiple::create(std::istream& istr) const {
  return std::make_shared<ComputeAddMultiple>(istr);
}

ComputeAddMultiple::ComputeAddMultiple(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeAddMultiple", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(834 == version, "mismatch version: " << version);
}

void ComputeAddMultiple::serialize_compute_add_multiple_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(834, ostr);
}

void ComputeAddMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_add_multiple_(ostr);
}

}  // namespace feasst
