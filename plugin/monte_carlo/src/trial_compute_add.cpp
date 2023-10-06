#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialComputeAdd::TrialComputeAdd(argtype args) : TrialComputeAdd(&args) {
  FEASST_CHECK_ALL_USED(args);
}
TrialComputeAdd::TrialComputeAdd(argtype * args) : TrialCompute(args) {
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
  const int iconf = stages->front()->select().configuration_index();
  acceptance->add_to_energy_new(criteria->current_energy(iconf), iconf);
  //acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile(iconf), iconf);
  acceptance->add_to_macrostate_shift(1, iconf);
  DEBUG("iconf " << iconf);
  DEBUG("old energy " << criteria->current_energy(iconf));
  DEBUG("new energy " << acceptance->energy_new(iconf));
  DEBUG("old energy profile " << feasst_str(criteria->current_energy_profile(iconf)));
  DEBUG("new energy profile " << feasst_str(acceptance->energy_profile_new(iconf)));
  { // Metropolis
    const Configuration& config = system->configuration(iconf);
    const double volume = config.domain().volume();
    const TrialSelect& select = (*stages)[0]->trial_select();
    const int particle_index = select.mobile().particle_index(0);
    const int particle_type = config.select_particle(particle_index).type();
    acceptance->set_macrostate_shift_type(particle_type, iconf);
    const double prob = 1./(config.num_particles_of_type(particle_type) + 1);
    DEBUG("volume " << volume << " prob " << prob << " betamu " << system->thermo_params().beta_mu(particle_type));
    acceptance->add_to_ln_metropolis_prob(
      std::log(volume*prob)
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
