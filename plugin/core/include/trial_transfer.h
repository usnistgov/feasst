
#ifndef FEASST_CORE_TRIAL_TRANSFER_H_
#define FEASST_CORE_TRIAL_TRANSFER_H_

#include "core/include/trial.h"
#include "core/include/perturb_transfer.h"
#include "core/include/random.h"
#include "core/include/criteria.h"

namespace feasst {

/**
 */
class TrialTransfer : public Trial {
 public:
  /// set the type of particle
  int particle_type_ = 0;

  void attempt(Criteria* criteria, System * system) {
    perturb_.before_attempt();
    criteria->before_attempt(system);
    if (random_.uniform() < add_probability_) {
      attempt_to_add(criteria, system);
    } else {
      attempt_to_remove(criteria, system);
    }
  }

  void attempt_to_remove(Criteria* criteria, System * system) {
    Configuration * config = system->configuration(0);
    /// HWH: change this up for group/type/controlled by config
    double delta_energy = 0;
    config->select_random_particle_of_type(particle_type_);
    if (config->selection().is_empty()) {
      accept_criteria_.force_rejection = 1;
    } else {
      delta_energy = -system->energy_of_selection();
      const int num_mol_old = config->num_particles();
      perturb_.remove_selected_particle(system);
      DEBUG("delta_energy " << delta_energy);
      accept_criteria_.ln_metropolis_prob =
        log(double(num_mol_old)/config->domain().volume())
        - criteria->beta()*delta_energy
        - log(criteria->activity());
      //-criteria->beta()*de;
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.force_rejection = 0;
      accept_criteria_.system = system;
    }
    if (criteria->is_accepted(accept_criteria_)) {
      DEBUG("accept del " << delta_energy);
    } else {
      DEBUG("reject del " << delta_energy);
      perturb_.revert();
    }
  }

  void attempt_to_add(Criteria* criteria, System * system) {
    Configuration * config = system->configuration(0);
    Particle particle = config->particle_types().particle(particle_type_);
    const Position rand_in_box = config->domain().random_position(&random_);
    DEBUG("rand_in_box " << rand_in_box.str());
    particle.displace(rand_in_box);
    perturb_.add(particle, system);
    config->select_last_particle();
    const int num_mol_new = config->num_particles();
    const double delta_energy = system->energy_of_selection();
    DEBUG("delta_energy " << delta_energy);
    accept_criteria_.ln_metropolis_prob =
      log(config->domain().volume()/double(num_mol_new))
      - criteria->beta()*delta_energy
      + log(criteria->activity());
    accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
    accept_criteria_.force_rejection = 0;
    accept_criteria_.system = system;
    if (criteria->is_accepted(accept_criteria_)) {
      DEBUG("accepted add " << delta_energy);
    } else {
      DEBUG("reject add " << delta_energy);
      perturb_.revert();
    }
  }

  void set_add_probability(const double prob) {
    add_probability_ = prob;
  }

  virtual ~TrialTransfer() {}

 private:
  PerturbTransfer perturb_;
  double add_probability_ = 0.5;
  Random random_;
  AcceptanceCriteria accept_criteria_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSFER_H_
