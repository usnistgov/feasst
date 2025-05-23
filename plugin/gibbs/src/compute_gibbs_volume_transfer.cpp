#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_select.h"
#include "gibbs/include/compute_gibbs_volume_transfer.h"

namespace feasst {

ComputeGibbsVolumeTransfer::ComputeGibbsVolumeTransfer() {
  class_name_ = "ComputeGibbsVolumeTransfer";
}

FEASST_MAPPER(ComputeGibbsVolumeTransfer,);

void ComputeGibbsVolumeTransfer::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeGibbsVolumeTransfer");
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  int config_del = 0;
  int config_add = 1;
  if ((*stages)[0]->trial_select().configuration_index() == 0) {
    config_add = 0;
    config_del = 1;
  }
  DEBUG("config_add " << config_add);
  DEBUG("config_del " << config_del);
  const Configuration& conf_add = system->configuration(config_add);
  const Configuration& conf_del = system->configuration(config_del);
  const double volume_old_add = conf_add.domain().volume();
  const double volume_old_del = conf_del.domain().volume();
  DEBUG("volume_old_add " << volume_old_add);
  DEBUG("volume_old_del " << volume_old_del);
  compute_rosenbluth(0, criteria, system, acceptance, stages, random);
  if (!acceptance->reject()) {
    for (int iconf : {config_add, config_del}) {
      acceptance->set_energy_new(acceptance->energy_new(iconf), iconf);
      acceptance->set_energy_profile_new(acceptance->energy_profile_new(iconf), iconf);
    }
    DEBUG("en 0 current " << MAX_PRECISION << criteria->current_energy(0));
    DEBUG("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
    DEBUG("en 1 current " << MAX_PRECISION << criteria->current_energy(1));
    DEBUG("en 1 new " << MAX_PRECISION << acceptance->energy_new(1));
    const double volume_new_add = conf_add.domain().volume();
    const double volume_new_del = conf_del.domain().volume();
    DEBUG("volume_new_add " << volume_new_add);
    DEBUG("volume_new_del " << volume_new_del);
    if (volume_old_add == volume_new_add) acceptance->set_reject(true);
    if (volume_old_del == volume_new_del) acceptance->set_reject(true);
    const int num_particles_from_add = conf_add.num_particles();
    const int num_particles_from_del = conf_del.num_particles();
    const ThermoParams& thermo = system->thermo_params();
    acceptance->add_to_ln_metropolis_prob(
      + num_particles_from_add*std::log(volume_new_add/volume_old_add)
      + num_particles_from_del*std::log(volume_new_del/volume_old_del)
      // manually add the energy of the old configurations
      // this is an optimization to avoid recomputing the old energy
      + thermo.beta()*criteria->current_energy(config_add)
      + thermo.beta()*criteria->current_energy(config_del)
    );
  }
}

std::shared_ptr<TrialCompute> ComputeGibbsVolumeTransfer::create(std::istream& istr) const {
  return std::make_shared<ComputeGibbsVolumeTransfer>(istr);
}

ComputeGibbsVolumeTransfer::ComputeGibbsVolumeTransfer(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeGibbsVolumeTransfer", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7345 == version, "mismatch version: " << version);
}

void ComputeGibbsVolumeTransfer::serialize_compute_gibbs_volume_transfer_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(7345, ostr);
}

void ComputeGibbsVolumeTransfer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_gibbs_volume_transfer_(ostr);
}

}  // namespace feasst
