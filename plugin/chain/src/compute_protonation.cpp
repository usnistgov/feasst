#include <math.h> // log10
#include "chain/include/compute_protonation.h"

namespace feasst {

ComputeProtonation::ComputeProtonation(const argtype& args) {
  class_name_ = "ComputeProtonation";
  Arguments args_(args);
  pKa_ = args_.key("pKa").dflt("0").dble();
}

class MapComputeProtonation {
 public:
  MapComputeProtonation() {
    auto obj = MakeComputeProtonation();
    obj->deserialize_map()["ComputeProtonation"] = obj;
  }
};

static MapComputeProtonation mapper_ = MapComputeProtonation();

void ComputeProtonation::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeProtonation");

  // Compute the old rosenbluth of both stages.
  compute_rosenbluth(1, criteria, system, acceptance, stages, random);

  // mid stage sets sites unphysical again
  for (TrialStage * stage : *stages) stage->mid_stage(system);

  // Finally, only compute the new rosenbluth of the first stage,
  // assuming that the second stage is when a particle is removed.
  if (stage0.size() == 0) {
    stage0.resize(1);
  }
  stage0[0] = (*stages)[0];
  DEBUG("computing new rosenbluth");
  compute_rosenbluth(0, criteria, system, acceptance, &stage0, random);

  // manually set the second stage physical again.
  (*stages)[1]->set_mobile_physical(true, system);

  const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
  acceptance->set_energy_new(criteria->current_energy() + delta_energy);
  DEBUG("delta_energy " << delta_energy);
  { // Metropolis
    const Configuration& config = system->configuration();
    const double volume = config.domain().volume();
    const TrialSelect * select_remove = (*stages)[1]->trial_select();
    DEBUG("to remove " << select_remove->mobile().str());
    const int particle_index = select_remove->mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    DEBUG("volume " << volume << " selprob " << select_remove->probability() << " betamu " << criteria->beta_mu(particle_type));
    double prob = volume*select_remove->probability();
    const TrialSelect * select_deprot = (*stages)[0]->trial_select();
    DEBUG("num_deprot " << round(1./select_deprot->probability()));
    DEBUG("select_deprot prob " << select_deprot->probability());
    prob /= select_deprot->probability();
    DEBUG("site index " << select_deprot->mobile().site_index(0, 0));
    const int prot_index = select_deprot->mobile().particle_index(0);
    const Particle& part = config.select_particle(prot_index);
    const int prot_type = part.site(select_deprot->mobile().site_index(0, 0)).type();
    DEBUG("prot_type " << prot_type);
    int num_prot = part.num_sites_of_type(prot_type);
    DEBUG("num_prot " << num_prot);
    if (num_prot == 0) {
      ASSERT(acceptance->reject(), "should have been rejected?");
      num_prot = 1;
    }
    prob /= static_cast<double>(num_prot);

    acceptance->add_to_ln_metropolis_prob(
      - log(prob)
      - criteria->beta_mu(particle_type)
      - log10(criteria->pH() - pKa_)
    );
  }
}

std::shared_ptr<TrialCompute> ComputeProtonation::create(std::istream& istr) const {
  return std::make_shared<ComputeProtonation>(istr);
}

ComputeProtonation::ComputeProtonation(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeProtonation", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(258 == version, "mismatch version: " << version);
  feasst_deserialize(&pKa_, istr);
}

void ComputeProtonation::serialize_compute_protonation_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(258, ostr);
  feasst_serialize(pKa_, ostr);
}

void ComputeProtonation::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_protonation_(ostr);
}

}  // namespace feasst
