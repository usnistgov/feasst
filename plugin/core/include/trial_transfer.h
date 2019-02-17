
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
  void attempt(Criteria* criteria, System * system) {
    if (random_.uniform() < add_probability_) {
      before_attempt(criteria, system, &add_);
      attempt_to_add(criteria, system);
    } else {
      before_attempt(criteria, system, &remove_);
      attempt_to_remove(criteria, system);
    }
  }

  void attempt_to_remove(Criteria* criteria, System * system) {
    Configuration * config = system->get_configuration();
    double delta_energy = 0;
    remove_.select_random_particle(particle_type_, config);
    if (remove_.selection().is_empty()) {
      accept_criteria_.force_rejection = 1;
    } else {
      delta_energy = -system->energy(remove_.selection());
      const int group = config->particle_type_to_group(particle_type_); // HWH optimize this
      const int num_mol_old = config->num_particles(group);
      remove_.remove_selected_particle(system);
      DEBUG("delta_energy " << delta_energy);
      accept_criteria_.ln_metropolis_prob =
        log(double(num_mol_old)/config->domain().volume())
        - criteria->beta()*delta_energy
        - log(criteria->activity());
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.force_rejection = 0;
      accept_criteria_.system = system;
    }
    accept_or_reject(accept_criteria_, &remove_, criteria);
  }

  void attempt_to_add(Criteria* criteria, System * system) {
    Configuration * config = system->get_configuration();
    const Position rand_in_box = config->domain().random_position(&random_);
    DEBUG("rand_in_box " << rand_in_box.str());
    add_.add(particle_type_, rand_in_box, system);
    const int num_mol_new = config->num_particles();
    const double delta_energy = system->energy(add_.selection());
    DEBUG("delta_energy " << delta_energy);
    accept_criteria_.ln_metropolis_prob =
      log(config->domain().volume()/double(num_mol_new))
      - criteria->beta()*delta_energy
      + log(criteria->activity());
    accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
    accept_criteria_.force_rejection = 0;
    accept_criteria_.system = system;
    accept_or_reject(accept_criteria_, &add_, criteria);
  }

  void set_add_probability(const double prob) {
    add_probability_ = prob;
  }

  virtual ~TrialTransfer() {}

 private:
  PerturbAdd add_;
  PerturbRemove remove_;
  double add_probability_ = 0.5;
  Random random_;
  AcceptanceCriteria accept_criteria_;

  /// set the type of particle
  int particle_type_ = 0;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSFER_H_
