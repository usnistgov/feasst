#include <math.h> // log10
#include "chain/include/compute_deprotonation.h"

namespace feasst {

ComputeDeprotonation::ComputeDeprotonation(const argtype& args) {
  class_name_ = "ComputeDeprotonation";
  Arguments args_(args);
  pKa_ = args_.key("pKa").dflt("0").dble();
}

class MapComputeDeprotonation {
 public:
  MapComputeDeprotonation() {
    auto obj = MakeComputeDeprotonation();
    obj->deserialize_map()["ComputeDeprotonation"] = obj;
  }
};

static MapComputeDeprotonation mapper_ = MapComputeDeprotonation();

void ComputeDeprotonation::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeDeprotonation");

  // To begin, only compute the old rosenbluth of the first stage,
  // assuming that the second stage is when a particle is added.
  if (stage0.size() == 0) {
    stage0.resize(1);
  }
  stage0[0] = (*stages)[0];
  compute_rosenbluth(1, criteria, system, acceptance, &stage0, random);

  // Compute the new rosenbluth of both stages.
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  const TrialSelect * select_add = (*stages)[1]->trial_select();
  if ((*stages)[0]->rosenbluth().num() > 1) {
    system->get_configuration()->revive(select_add->mobile());
  }
  const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();
    const int particle_index = select_add->mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    DEBUG("volume " << volume << " selprob " << select_add->probability() << " betamu " << criteria->beta_mu(particle_type));
    double prob = volume*select_add->probability();
    const TrialSelect * select_prot = (*stages)[0]->trial_select();
    DEBUG("num_prot " << round(1./select_prot->probability()));
    DEBUG("select_prot prob " << select_prot->probability());
    prob /= select_prot->probability();
    DEBUG("site index " << select_prot->mobile().site_index(0, 0));
    const int deprot_index = select_prot->mobile().particle_index(0);
    const Particle& part = config.select_particle(deprot_index);
    const int deprot_type = part.site(select_prot->mobile().site_index(0, 0)).type();
    DEBUG("deprot_type " << deprot_type);
    int num_deprot = part.num_sites_of_type(deprot_type);
    DEBUG("num_deprot " << num_deprot);
    if (num_deprot == 0) {
      ASSERT(acceptance->reject(), "should have been rejected?");
      num_deprot = 1;
    }
    prob /= static_cast<double>(num_deprot);

    acceptance->add_to_ln_metropolis_prob(
      log(prob)
      + criteria->beta_mu(particle_type)
      + log10(criteria->pH() - pKa_)
    );
  }
}

std::shared_ptr<TrialCompute> ComputeDeprotonation::create(std::istream& istr) const {
  return std::make_shared<ComputeDeprotonation>(istr);
}

ComputeDeprotonation::ComputeDeprotonation(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeDeprotonation", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(258 == version, "mismatch version: " << version);
  feasst_deserialize(&pKa_, istr);
}

void ComputeDeprotonation::serialize_compute_deprotonation_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(258, ostr);
  feasst_serialize(pKa_, ostr);
}

void ComputeDeprotonation::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_deprotonation_(ostr);
}

}  // namespace feasst
