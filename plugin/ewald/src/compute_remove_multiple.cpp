#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "ewald/include/compute_remove_multiple.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

ComputeRemoveMultiple::ComputeRemoveMultiple(const argtype& args) {
  class_name_ = "ComputeRemoveMultiple";
  Arguments args_(args);
  args_.dont_check();
  shift_ = args_.key("shift").dflt("-1").integer();
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
  acceptance->add_to_macrostate_shift(shift_);
  acceptance->set_macrostate_shift_type(-1);  // disable constraints with multi-particles
  DEBUG("deltaE " << MAX_PRECISION << -1.*acceptance->energy_old());
  const Configuration& config = system->configuration();
  const double volume = config.domain().volume();

  // initialize delta_
  const int num_ptypes = config.num_particle_types();
  if (static_cast<int>(delta_.size()) < num_ptypes) delta_.resize(num_ptypes);
  std::fill(delta_.begin(), delta_.end(), 0);

  for (const TrialStage * stage : *stages) {
    const TrialSelect& select = stage->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    DEBUG("volume " << volume << " selprob " << select.probability()
      << " betamu " << criteria->beta_mu(particle_type));
    const int num_pt = config.num_particles_of_type(particle_type);
    const double prob = 1./static_cast<double>(num_pt - delta_[particle_type]);
    ++delta_[particle_type];
    acceptance->add_to_ln_metropolis_prob(
      - std::log(volume*prob)
      //- std::log(volume*select.probability())
      - criteria->beta_mu(particle_type)
    );
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
  feasst_deserialize(&shift_, istr);
}

void ComputeRemoveMultiple::serialize_compute_remove_multiple_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(294, ostr);
  feasst_serialize(shift_, ostr);
}

void ComputeRemoveMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_remove_multiple_(ostr);
}

}  // namespace feasst
