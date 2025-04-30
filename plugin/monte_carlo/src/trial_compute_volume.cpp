#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialComputeVolume::TrialComputeVolume(argtype args) : TrialCompute(&args) {
  class_name_ = "TrialComputeVolume";
  feasst_check_all_used(args);
}

FEASST_MAPPER(TrialComputeVolume,);

void TrialComputeVolume::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeVolume");
  const int iconf = (*stages)[0]->select().configuration_index();
  DEBUG("iconf " << iconf);
  const Configuration& config = system->configuration(iconf);
  const double volume_old = config.domain().volume();
  DEBUG("volume_old " << volume_old);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  if (!acceptance->reject()) {
    acceptance->set_energy_new(acceptance->energy_new(iconf), iconf);
    DEBUG("energy new " << acceptance->energy_new(iconf));
    acceptance->set_energy_profile_new(acceptance->energy_profile_new(iconf), iconf);
    DEBUG("acceptance en prof " << feasst_str(acceptance->energy_profile_new(iconf)));
    const double volume_new = config.domain().volume();
    DEBUG("volume_new " << volume_new);
    if (volume_old == volume_new) acceptance->set_reject(true);
    const ThermoParams& thermo = system->thermo_params();
    DEBUG("dv " << volume_new - volume_old);
    acceptance->add_to_ln_metropolis_prob(
      - thermo.beta()*thermo.pressure()*(volume_new - volume_old)
      + config.num_particles()*std::log(volume_new/volume_old)
      // manually add the energy of the old configuration
      // this is an optimization to avoid recomputing the old energy
      + thermo.beta()*criteria->current_energy(iconf)
    );
  }
}

std::shared_ptr<TrialCompute> TrialComputeVolume::create(std::istream& istr) const {
  return std::make_shared<TrialComputeVolume>(istr);
}

TrialComputeVolume::TrialComputeVolume(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "TrialComputeVolume", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(834 == version, "mismatch version: " << version);
}

void TrialComputeVolume::serialize_trial_compute_add_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(834, ostr);
}

void TrialComputeVolume::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_add_(ostr);
}

}  // namespace feasst
