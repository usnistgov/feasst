
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
  TrialTransfer(const argtype &args = argtype()) : Trial(args) {}

  void attempt(Criteria* criteria, System * system) override {
    if (random_.uniform() < add_probability_) {
      attempt_to_add(criteria, system);
    } else {
      attempt_to_remove(criteria, system);
    }
    num_particles_of_type_ =
      system->configuration().num_particles_of_type(particle_type_);
    DEBUG("finished " << add_probability_);
  }

  void attempt_to_remove(Criteria* criteria, System * system) {
    DEBUG("attempt to remove");
    before_attempt(criteria, system, &remove_, &accept_criteria_);
    Configuration * config = system->get_configuration();
    double delta_energy = 0;
    remove_.select_random_particle(particle_type_, config);
    if (remove_.selection().is_empty()) {
      accept_criteria_.force_rejection = 1;
    } else {
      delta_energy = -system->energy(remove_.selection());
      const int num_mol_old = config->num_particles_of_type(particle_type_);
      remove_.remove_selected_particle(system);
      DEBUG("delta_energy " << delta_energy);
      accept_criteria_.ln_metropolis_prob =
        log(double(num_mol_old)/config->domain().volume())
        - criteria->beta()*delta_energy
        - criteria->beta_mu(particle_type_);
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.force_rejection = 0;
      accept_criteria_.system = system;
    }
    accept_or_reject(accept_criteria_, &remove_, criteria);
  }

  void attempt_to_add(Criteria* criteria, System * system) {
    DEBUG("attempt to add");
    before_attempt(criteria, system, &add_, &accept_criteria_);
    Configuration * config = system->get_configuration();
    const Position rand_in_box = config->domain().random_position(&random_);
    DEBUG("rand_in_box " << rand_in_box.str());
    add_.add(particle_type_, rand_in_box, system);
    const int num_mol_new = config->num_particles(particle_type_);
    const double delta_energy = system->energy(add_.selection());
    DEBUG("delta_energy " << delta_energy);
    accept_criteria_.ln_metropolis_prob =
      log(config->domain().volume()/double(num_mol_new))
      - criteria->beta()*delta_energy
      + criteria->beta_mu(particle_type_);
    accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
    accept_criteria_.force_rejection = 0;
    accept_criteria_.system = system;
    accept_or_reject(accept_criteria_, &add_, criteria);
  }

  void set_add_probability(const double prob) {
    add_probability_ = prob;
  }

  // see base class
  std::string status_header() const override {
    std::stringstream ss;
    ss << Trial::status_header() << " N" << particle_type_;
    return ss.str();
  }

  // see base class
  std::string status() const override {
    std::stringstream ss;
    ss << Trial::status() << " " << num_particles_of_type_;
    return ss.str();
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

  /// store the number of particles of type
  int num_particles_of_type_ = 0;
};

inline std::shared_ptr<TrialTransfer> MakeTrialTransfer(const argtype &args = argtype()) {
  return std::make_shared<TrialTransfer>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSFER_H_
