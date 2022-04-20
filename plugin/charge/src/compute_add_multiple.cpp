#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select.h"
#include "charge/include/compute_add_multiple.h"

namespace feasst {

ComputeAddMultiple::ComputeAddMultiple(argtype args)
  : ComputeAddMultiple(&args) {
  check_all_used(args);
}
ComputeAddMultiple::ComputeAddMultiple(argtype * args)
  : TrialCompute(args) {
  class_name_ = "ComputeAddMultiple";
  shift_ = integer("shift", args, 1);
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
  acceptance->add_to_energy_new(criteria->current_energy());
  //acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile());
  acceptance->add_to_macrostate_shift(shift_);
  acceptance->set_macrostate_shift_type(-1);  // disable constraints with multi-particles
  DEBUG("deltaE " << MAX_PRECISION << acceptance->energy_new());
  const Configuration& config = system->configuration();
  const double volume = config.domain().volume();

  // initialize delta_
  const int num_ptypes = config.num_particle_types();
  if (static_cast<int>(delta_.size()) < num_ptypes) delta_.resize(num_ptypes);
  std::fill(delta_.begin(), delta_.end(), 0);

  DEBUG(config.num_particles());
  DEBUG(feasst_str(delta_));

  // Metropolis
  for (const TrialStage * stage : *stages) {
    DEBUG("stage");
    const TrialSelect& select = stage->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    const int num_pt = config.num_particles_of_type(particle_type);
    ++delta_[particle_type];
    const double prob = 1./static_cast<double>(num_pt + delta_[particle_type]);
    DEBUG("volume " << volume << " prob " << prob <<
      " betamu " << system->thermo_params().beta_mu(particle_type));
    acceptance->add_to_ln_metropolis_prob(
      std::log(volume*prob)
      + system->thermo_params().beta_mu(particle_type)
    );
  }
}

std::shared_ptr<TrialCompute> ComputeAddMultiple::create(std::istream& istr) const {
  return std::make_shared<ComputeAddMultiple>(istr);
}

ComputeAddMultiple::ComputeAddMultiple(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeAddMultiple", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(9346 == version, "mismatch version: " << version);
  feasst_deserialize(&shift_, istr);
}

void ComputeAddMultiple::serialize_compute_add_multiple_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(9346, ostr);
  feasst_serialize(shift_, ostr);
}

void ComputeAddMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_add_multiple_(ostr);
}

}  // namespace feasst
