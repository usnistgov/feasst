#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select.h"
#include "gibbs/include/compute_gibbs_particle_transfer.h"

namespace feasst {

ComputeGibbsParticleTransfer::ComputeGibbsParticleTransfer() {
  class_name_ = "ComputeGibbsParticleTransfer";
}

class MapComputeGibbsParticleTransfer {
 public:
  MapComputeGibbsParticleTransfer() {
    auto obj = MakeComputeGibbsParticleTransfer();
    obj->deserialize_map()["ComputeGibbsParticleTransfer"] = obj;
  }
};

static MapComputeGibbsParticleTransfer mapper_ = MapComputeGibbsParticleTransfer();

void ComputeGibbsParticleTransfer::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("ComputeGibbsParticleTransfer");
  DEBUG("lnmet " << acceptance->ln_metropolis_prob());

  // del
  std::vector<TrialStage*> del_stages = {(*stages)[1]};
  compute_rosenbluth(1, criteria, system, acceptance, &del_stages, random);
  std::vector<TrialStage*> add_stages = {(*stages)[0]};
  compute_rosenbluth(0, criteria, system, acceptance, &add_stages, random);
  int config_del = 0;
  int config_add = 1;
  if ((*stages)[0]->trial_select().configuration_index() == 0) {
    config_del = 1;
    config_add = 0;
  }
  DEBUG("config_add " << config_add);
  DEBUG("config_del " << config_del);
  acceptance->set_energy_new(criteria->current_energy(config_del) - acceptance->energy_old(config_del), config_del);
  acceptance->set_energy_profile_new(criteria->current_energy_profile(config_del), config_del);
  acceptance->subtract_from_energy_profile_new(acceptance->energy_profile_old(config_del), config_del);
  acceptance->add_to_macrostate_shift(-1, config_del);

  DEBUG("energy contribution of config " << config_del << " particle to be deleted: " << acceptance->energy_old(config_del));

  // add
  acceptance->add_to_energy_new(criteria->current_energy(config_add), config_add);
  acceptance->add_to_energy_profile_new(criteria->current_energy_profile(config_add), config_add);
  acceptance->add_to_macrostate_shift(1, config_add);

  DEBUG("energy contribution of config " << config_add << " particle to be added: " << acceptance->energy_new(config_add));

  DEBUG("en 0 current " << MAX_PRECISION << criteria->current_energy(0));
  DEBUG("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
  DEBUG("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));
  DEBUG("en 1 current " << MAX_PRECISION << criteria->current_energy(1));
  DEBUG("en 1 new " << MAX_PRECISION << acceptance->energy_new(1));
  //DEBUG("en 1 old acc " << MAX_PRECISION << acceptance->energy_old(1));

  DEBUG("lnmet " << acceptance->ln_metropolis_prob());
//  DEBUG("en 0 old " << criteria->current_energy(0));
//  DEBUG("en 0 new " << MAX_PRECISION << acceptance->energy_new(0));
//  DEBUG("en 0 old acc " << MAX_PRECISION << acceptance->energy_old(0));
//  DEBUG("en 1 old " << criteria->current_energy(1));
//  DEBUG("en 1 new " << MAX_PRECISION << acceptance->energy_new(1));
//  DEBUG("en 1 old acc " << MAX_PRECISION << acceptance->energy_old(1));
  { // Metropolis
    const Configuration& conf_add = system->configuration(config_add);
    const Configuration& conf_del = system->configuration(config_del);
    const TrialSelect& select_add = (*stages)[0]->select();
    //const TrialSelect& select_del = (*stages)[1]->select();
    const int particle_add = select_add.mobile().particle_index(0);
    //const int particle_del = select_del.mobile().particle_index(0);
    const int particle_type = conf_add.select_particle(particle_add).type();
    acceptance->set_macrostate_shift_type(particle_type, config_add);
    acceptance->set_macrostate_shift_type(particle_type, config_del);
//    DEBUG("lnselprob " << std::log(select.probability()));
//    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    const int num_particles_from_add = conf_add.num_particles_of_type(particle_type);
    const int num_particles_from_del = conf_del.num_particles_of_type(particle_type);
    DEBUG("num_particles_from_add " << num_particles_from_add);
    DEBUG("num_particles_from_del " << num_particles_from_del);
    const double vol_from_add = conf_add.domain().volume();
    const double vol_from_del = conf_del.domain().volume();
    DEBUG("vol_from_add " << vol_from_add);
    DEBUG("vol_from_del " << vol_from_del);
    acceptance->add_to_ln_metropolis_prob(
      std::log(num_particles_from_del/static_cast<double>(num_particles_from_add+1)*
        vol_from_add/vol_from_del)
      //std::log(select.probability())
    );
    DEBUG("lnmet " << acceptance->ln_metropolis_prob());
  }
}

std::shared_ptr<TrialCompute> ComputeGibbsParticleTransfer::create(std::istream& istr) const {
  return std::make_shared<ComputeGibbsParticleTransfer>(istr);
}

ComputeGibbsParticleTransfer::ComputeGibbsParticleTransfer(std::istream& istr)
  : TrialCompute(istr) {
  // ASSERT(class_name_ == "ComputeGibbsParticleTransfer", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7345 == version, "mismatch version: " << version);
}

void ComputeGibbsParticleTransfer::serialize_compute_gibbs_particle_transfer_(std::ostream& ostr) const {
  serialize_trial_compute_(ostr);
  feasst_serialize_version(7345, ostr);
}

void ComputeGibbsParticleTransfer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_gibbs_particle_transfer_(ostr);
}

}  // namespace feasst
