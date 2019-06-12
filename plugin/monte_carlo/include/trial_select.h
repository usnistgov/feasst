
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_H_

#include <vector>
#include <numeric>
#include <string>
#include <memory>
#include "system/include/system.h"
#include "system/include/select_list.h"
#include "monte_carlo/include/criteria.h"
#include "utils/include/arguments.h"

namespace feasst {

/**
  Select the mobile particles and sites that are to be perturbed via trials.
  Store the original position in mobile_original for reverting.
  Also store the 'anchor' particles and sites which are not mobile but may be
  required to complete the perturbation (e.g., bonds).
 */
class TrialSelect {
 public:
  TrialSelect(
    /**
      group_index : index of group definied within system (default: 0).

      particle_type : type of particle in configuration (default: -1)
     */
    const argtype& args = argtype()) {
    // set_mayer();
    args_.init(args);
    args_.dont_check();
    group_index_ = args_.key("group_index").dflt("0").integer();
    particle_type_ = args_.key("particle_type").dflt("0").integer();
    set_probability();
  }

  /// Return the index of group for selection.
  int group_index() const { return group_index_; }

  /// Return the particle type.
  int particle_type() const {
    ASSERT(particle_type_ != -1, "particle type not specified");
    return particle_type_;
  }

  /// Perform the selection as implemented in the derived class.
  virtual void select(System * system) { ERROR("not implemented"); }

  /// Check if the selection was completed properly, or if the trial should
  /// be automatically rejected.
  virtual void check(Acceptance * acceptance) const {
    if (mobile_.num_sites() == 0) {
      acceptance->set_reject(true);
    }
  }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(System * system) {};

  /// Return the mobile selection. These can change during the trial.
  const SelectList& mobile() const { return mobile_; }

  /// Return a pointer to the mobile selection.
  SelectList * get_mobile() { return &mobile_; }

  /// Set the mobile selection.
  void set_mobile(const SelectList& mobile) { mobile_ = mobile; }

  /// Return originally-seleted mobile. These do not change during trial.
  const SelectList& mobile_original() const { return mobile_original_; }

  /// Return the anchor selection.
  const SelectList& anchor() const { return anchor_; }

  /// Set the state of the trial for the mobile select (e.g., old, move, add").
  void set_trial_state(const std::string state) {
    mobile_.set_trial_state(state); }

  /// Reset the mobile selection to the original.
  void reset_mobile() { mobile_ = mobile_original_; }

  /// Return the probability of the selection. For example, if a random particle
  /// type is selected, then the probability is the inverse of the number of
  /// particles of that type.
  double probability() const { return probability_; }

  /// Set the probability of selection.
  void set_probability(const double prob = 1) { probability_ = prob; }

 protected:
  SelectList mobile_original_;
  SelectList mobile_;
  SelectList anchor_;

 private:
  Arguments args_;
  int group_index_;
  int particle_type_;

  // optimzation or temporary object
  double probability_;
};

/// Select a random particle for trial.
class TrialSelectParticle : public TrialSelect {
 public:
  TrialSelectParticle(
    /**
      load_coordinates : load the coordinates into the selection (default: true)
     */
    const argtype& args = argtype()) : TrialSelect(args) {
    Arguments args_(args);
    args_.dont_check();
    load_coordinates_ = args_.key("load_coordinates").dflt("true").boolean();
  }

  bool load_coordinates() const { return load_coordinates_; }

  void select(System* system) override {
    // HWH consider removing this from SelectList and putting here instead.
    int load = 0;
    if (load_coordinates()) load = 1;
    const int num = mobile_.random_particle(system->configuration(),
      group_index(),
      load);
    set_probability(1./static_cast<double>(num));
    mobile_original_ = mobile_;
  }

 private:
  bool load_coordinates_;
};

inline std::shared_ptr<TrialSelectParticle> MakeTrialSelectParticle(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectParticle>(args);
}

/// Select a random particle of a given type.
class TrialSelectParticleOfType : public TrialSelectParticle {
 public:
  TrialSelectParticleOfType(const argtype& args = argtype())
    : TrialSelectParticle(args) {}
  void select(System* system) override {
    // HWH consider removing this from SelectList and putting here instead.
    int load = 0;
    if (load_coordinates()) load = 1;
    const int num = mobile_.random_particle_of_type(
      particle_type(),
      system->get_configuration(),
      load);
    set_probability(1./static_cast<double>(num));
    mobile_original_ = mobile_;
  }
};

/**
  The selection in a trial which adds particles, for example, does nothing.
  Rather, the perturbation adds the particle and selects the added one.
 */
class TrialSelectDoNothing : public TrialSelect {
 public:
  void select(System* system) override {}
  void check(Acceptance * acceptance) const override {
    // HWH removed, because selection doesn't reset each trial for optimization
    //ASSERT(mobile_.num_sites() == 0, "Trials which add should not have a "
    //  << "selection at this stage.");
  }
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_H_
