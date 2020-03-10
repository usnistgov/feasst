#include "monte_carlo/include/trial_select.h"
#include "ewald/include/compute_remove_multiple.h"

namespace feasst {

ComputeRemoveMultiple::ComputeRemoveMultiple() {
  class_name_ = "ComputeRemoveMultiple";
}

class MapComputeRemoveMultiple {
 public:
  MapComputeRemoveMultiple() {
    auto obj = MakeComputeRemoveMultiple();
    obj->deserialize_map()["ComputeRemoveMultiple"] = obj;
  }
};

static MapComputeRemoveMultiple mapper_ = MapComputeRemoveMultiple();

void ComputeRemoveMultiple::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeRemoveMultiple");
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);
  acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
  acceptance->add_to_macrostate_shift(-1);
  acceptance->set_macrostate_shift_type(-1);  // disable constraints with multi-particles
  DEBUG("deltaE " << MAX_PRECISION << -1.*acceptance->energy_old());
  for (const TrialStage * stage : *stages) {
    { // Metropolis
      const Configuration& config = system->configuration();
      const double volume = config.domain()->volume();
      const TrialSelect * select = stage->trial_select();
      const int particle_index = select->mobile().particle_index(0);
      const int particle_type = config.select_particle(particle_index).type();
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
      DEBUG("volume " << volume << " selprob " << select->probability()
        << " betamu " << criteria->beta_mu(particle_type));
      acceptance->add_to_ln_metropolis_prob(
        - log(volume*select->probability())
        - criteria->beta_mu(particle_type)
      );
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    }
  }
}

std::shared_ptr<TrialCompute> ComputeRemoveMultiple::create(
    std::istream& istr) const {
  return std::make_shared<ComputeRemoveMultiple>(istr);
}

ComputeRemoveMultiple::ComputeRemoveMultiple(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeRemoveMultiple", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(294 == version, "mismatch version: " << version);
}

void ComputeRemoveMultiple::serialize_compute_remove_multiple_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(294, ostr);
}

void ComputeRemoveMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_remove_multiple_(ostr);
}

}  // namespace feasst
