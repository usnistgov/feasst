
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
class TrialSelect : public PropertiedEntity {
 public:
  TrialSelect(
    /**
      group_index : index of group definied within system (default: 0).

      particle_type : type of particle in configuration (default: -1)
     */
    const argtype& args = argtype()) : PropertiedEntity() {
    // set_mayer();
    args_.init(args);
    args_.dont_check();

    // parse particle type and group index from args.
    particle_type_ = 0;
    group_index_ = 0;
    if (args_.key("particle_type").used()) {
      particle_type_ = args_.integer();
      ASSERT(args_.key("group_index").empty(),
        "cant specify both particle type and group index");
    } else {
      if (args_.key("group_index").used()) {
        group_index_ = args_.key("group_index").integer();
      }
    }

    set_probability();
  }

  /// Return the index of group for selection.
  int group_index() const { return group_index_; }

  /// Return the particle type.
  int particle_type() const {
    ASSERT(particle_type_ != -1, "particle type not specified");
    return particle_type_;
  }

  /// Perform upkeep before select.
  void before_select() {
    mobile_.reset_excluded_and_bond();
  }

  /// Perform the selection as implemented in the derived class.
  virtual void select(System * system) { ERROR("not implemented"); }

  /// Check if the selection was completed properly, or if the trial should
  /// be automatically rejected.
  // HWH depreciate check, make check return if failed
  virtual void check(Acceptance * acceptance) const {
    if (mobile_.num_sites() == 0) {
      acceptance->set_reject(true);
    }
  }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(System * system) {}

  /// Return the mobile selection. These can change during the trial.
  const SelectList& mobile() const { return mobile_; }

  /// Return a pointer to the mobile selection.
  SelectList * get_mobile() { return &mobile_; }

  /// Set the mobile selection.
  void set_mobile(const SelectList& mobile) { mobile_ = mobile; }

  /// Return originally-seleted mobile. These do not change during trial.
  const SelectList& mobile_original() const { return mobile_original_; }

  /// Return the anchor selection.
  const Select& anchor() const { return anchor_; }

  /// Return anchor position.
  const Position& anchor_position(
    /// anchor index, not configuration index
    const int particle_index,
    /// anchor index
    const int site_index,
    const System * system) {
    const int part = anchor_.particle_index(particle_index);
    const int site = anchor_.site_index(particle_index, site_index);
    DEBUG("site " << site);
    return system->configuration().select_particle(part).site(site).position();
  }

  /// Set the state of the trial for the mobile select (e.g., old, move, add).
  /// See Select::trial_state
  void set_trial_state(const int state) {
    mobile_.set_trial_state(state); }

  /// Reset the mobile selection to the original.
  void reset_mobile() { mobile_ = mobile_original_; }

  /// Return the probability of the selection. For example, if a random particle
  /// type is selected, then the probability is the inverse of the number of
  /// particles of that type.
  double probability() const { return probability_; }

  /// Set the probability of selection.
  void set_probability(const double prob = 1) { probability_ = prob; }

  /// Call after old configuration but before new.
  virtual void mid_stage() {}

 protected:
  SelectList mobile_original_;
  SelectList mobile_;
  Select anchor_;

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

// HWH remove "Of" in name
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

/// Select a random particle of a given type, and then a particular site.
class TrialSelectSiteInParticleType : public TrialSelect {
 public:
  TrialSelectSiteInParticleType(
    /**
      site : site index to select (default: 0).
     */
    const argtype& args = argtype())
    : TrialSelect(args) {
    Arguments args_(args);
    args_.dont_check();
    const int site = args_.key("site").dflt("0").integer();
    mobile_.clear();
    mobile_.add_site(0, site);
  }

  virtual void precompute(System * system) {
    group_index_ = system->get_configuration()->particle_type_to_group(particle_type());
  }

  void select(System* system) override {
    const Configuration& config = system->configuration();
    const int num = config.num_particles(group_index_);
    if (num > 0) {
      const int particle = random_.uniform(0, num - 1);
      mobile_.set_particle(0, particle);
      mobile_.load_positions(config.particles());
      mobile_original_ = mobile_;
    }
  }

 private:
  int group_index_;
  Random random_;
};

// requires a particle type and the indices of the two sites (mobile and anchor)
// from the two sites, we can get a bond index.
// from the bond index, we can get bond type (for bond parameters, etc)
// but how do we pass bond length, etc, to the perturbation?
class TrialSelectBond : public TrialSelect {
 public:
  TrialSelectBond(
    /**
      mobile_site : index of the mobile site.
      anchor_site : index of the anchor site.
     */
    const argtype& args = argtype()) : TrialSelect(args) {
    Arguments args_(args);
    args_.dont_check();
    mobile_site_ = args_.key("mobile_site").integer();
    anchor_site_ = args_.key("anchor_site").integer();
  }

  void precompute(System * system) override {
    const Particle& part = system->configuration().particle_types().particle(particle_type());
    const int bond_type = part.bond(mobile_site_, anchor_site_).type();
    const Bond& bond = system->configuration().unique_types().particle(particle_type()).bond(bond_type);
    add_property("bond_length", bond.property("length"));
    anchor_.clear();
    anchor_.add_site(0, anchor_site_);
    mobile_.clear();
    mobile_.add_site(0, mobile_site_);
  }

  void select(System * system) override {
    // select random particle of correct type
    Configuration * config = system->get_configuration();
    const int group_index = config->particle_type_to_group(particle_type());
    const int num = config->num_particles(group_index);
    ASSERT(num > 0, "code in selection failure without clearning mobile");
    const int index = random_.uniform(0, num - 1);
    const SelectGroup& select = config->group_select(group_index);
    mobile_.set_particle(0, select.particle_index(index));
    mobile_.load_positions(config->particles());
    mobile_original_ = mobile_;
  }

 private:
  int mobile_site_;
  int anchor_site_;
  Random random_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_H_
