#include "chain/include/compute_synthesis.h"

namespace feasst {

ComputeSynthesis::ComputeSynthesis() {
  class_name_ = "ComputeSynthesis";
}

class MapComputeSynthesis {
 public:
  MapComputeSynthesis() {
    auto obj = MakeComputeSynthesis();
    obj->deserialize_map()["ComputeSynthesis"] = obj;
  }
};

static MapComputeSynthesis mapper_ = MapComputeSynthesis();

void ComputeSynthesis::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeSynthesis");
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  const TrialSelect * select = (*stages)[1]->trial_select();
  if ((*stages)[0]->rosenbluth().num() > 1) {
    system->get_configuration()->revive(select->mobile());
  }
  acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();
    const int particle_index = select->mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    DEBUG("volume " << volume << " selprob " << select->probability() << " betamu " << criteria->beta_mu(particle_type));
    acceptance->add_to_ln_metropolis_prob(
      log(volume*select->probability())
      + criteria->beta_mu(particle_type)
    );
  }
}

std::shared_ptr<TrialCompute> ComputeSynthesis::create(std::istream& istr) const {
  return std::make_shared<ComputeSynthesis>(istr);
}

ComputeSynthesis::ComputeSynthesis(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeSynthesis", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(258 == version, "mismatch version: " << version);
}

void ComputeSynthesis::serialize_compute_synthesis_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(258, ostr);
}

void ComputeSynthesis::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_synthesis_(ostr);
}

}  // namespace feasst
