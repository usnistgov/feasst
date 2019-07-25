
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
    particle_type_ = -1;
    group_index_ = 0;
    if (args_.key("particle_type").used()) {
      is_particle_type_set_ = true;
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
    ASSERT(is_particle_type_set_, "particle type not specified");
    return particle_type_;
  }

  /// Perform upkeep before select.
  void before_select() {
    mobile_.reset_excluded_and_bond();
  }

  /// Perform the selection as implemented in the derived class.
  /// Return false if the selection cannot be made. Otherwise, return true.
  virtual bool select(
    /// Perturbed is included to allow chaining of selection based on previous.
    const Select& perturbed,
    System * system) {
    ERROR("not implemented"); }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(System * system) {
    if (is_particle_type_set_) {
      group_index_ = system->get_configuration()->particle_type_to_group(
        particle_type_);
    }
  }

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

  /// Select from ghost particles.
  void set_ghost(const bool ghost = true) {
    is_ghost_ = ghost;
    if (is_ghost_) {
      ASSERT(group_index() == 0, "ghost particles cannot be selected by groups");
      ASSERT(particle_type() >= 0, "ghost particles must be selected by type");
    }
  }

  /// Return true if selecting from ghost particles.
  bool is_ghost() const { return is_ghost_; }

 protected:
  SelectList mobile_original_;
  SelectList mobile_;
  Select anchor_;
  Random * random() { return &random_; }

 private:
  Arguments args_;
  int group_index_;
  int particle_type_;
  bool is_particle_type_set_ = false;
  bool is_ghost_;

  // optimzation or temporary object
  double probability_;
  Random random_;
};

/// Select a random particle for trial.
class TrialSelectParticle : public TrialSelect {
 public:
  TrialSelectParticle(
    /**
      load_coordinates : load the coordinates into the selection (default: true)

      site : site index to select. If all sites, set to -1 (default).

      ghost : select ghost particles (default: false).
     */
    const argtype& args = argtype()) : TrialSelect(args) {
    Arguments args_(args);
    args_.dont_check();
    load_coordinates_ = args_.key("load_coordinates").dflt("true").boolean();

    // parse site
    site_ = args_.key("site").dflt("-1").integer();
    if (site_ != -1) {
      mobile_.clear();
      mobile_.add_site(0, site_);
      site_vec_ =  {site_};
    }

    set_ghost(args_.key("ghost").dflt("false").boolean());
  }

  /// Return true if loading coordinates into selection.
  bool load_coordinates() const { return load_coordinates_; }

  /// Add random particle in group index to select.
  /// Return the number of particles to choose from.
  int random_particle(const Configuration& config,
      SelectPosition * select) {
    ASSERT(group_index() >= 0, "error");
    const int num = config.num_particles(group_index());
    if (num > 0) {
      const int index = random()->uniform(0, num - 1);
      const SelectGroup& ran = config.group_select(group_index());
      DEBUG("index " << group_index() << " " << index);
      DEBUG("num " << ran.num_particles());
      bool fast;
      if (site_ == - 1) {
        fast = select->replace_indices(ran.particle_index(index),
                                       ran.site_indices(index));
      } else {
        fast = select->replace_indices(ran.particle_index(index),
                                       site_vec_);
      }
      if (load_coordinates()) {
        if (!fast) select->resize();
        select->load_positions(config.particles());
      }
    } else {
      select->clear();
    }
    return num;
  }

  /// Add random particle in group index to select.
  /// Return the number of particles to choose from.
  void ghost_particle(Configuration * config,
    SelectPosition * select) {
    ASSERT(static_cast<int>(config->ghosts().size()) > particle_type(),
      "type not recognized");
    // if no ghosts, create one
    DEBUG("num ghosts " << config->ghosts()[particle_type()].num_particles());
    if (config->ghosts()[particle_type()].num_particles() == 0) {
      config->add_particle_of_type(particle_type());
      Select add;
      add.add_particle(config->newest_particle(), config->newest_particle_index());
      config->remove_particles(add);
    }
    const Select& ghost = config->ghosts()[particle_type()];
    bool fast;
    // replace indices with the last ghost as may be optimal method available
    // to delete.
    if (site_ == -1) {
      fast = select->replace_indices(ghost.particle_indices().back(),
                                     ghost.site_indices().back());
    } else {
      fast = select->replace_indices(ghost.particle_indices().back(),
                                     site_vec_);
      config->set_selection_physical(ghost, false);
      config->set_selection_physical(*select, true);
    }
    if (load_coordinates()) {
      if (!fast) select->resize();
      select->load_positions(config->particles());
    }
  }

  bool select(const Select& perturbed, System* system) override {
    if (is_ghost()) {
      ghost_particle(system->get_configuration(), &mobile_);
      set_probability(1.);
    } else {
      const int num = random_particle(system->configuration(), &mobile_);
      if (num <= 0) return false;
      set_probability(1./static_cast<double>(num));
    }
    mobile_.remove_unphysical_sites(system->configuration());
    mobile_original_ = mobile_;
    return true;
  }

 private:
  bool load_coordinates_;
  int site_;
  std::vector<int> site_vec_; // optimization
};

inline std::shared_ptr<TrialSelectParticle> MakeTrialSelectParticle(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectParticle>(args);
}

/**
  A random particle of given type is selected if previously perturbed sites
    are not available.
  Select a single bond from given anchor to mobile sites.
 */
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

  int anchor_site() const { return anchor_site_; }
  int mobile_site() const { return mobile_site_; }

  // bond_length is added as a property
  // mobile and anchor are sized
  void precompute(System * system) override {
    TrialSelect::precompute(system);
    const Particle& part = system->configuration().particle_types().particle(particle_type());
    const int bond_type = part.bond(mobile_site_, anchor_site_).type();
    const Bond& bond = system->configuration().unique_types().particle(particle_type()).bond(bond_type);
    add_property("bond_length", bond.property("length"));
    anchor_.clear();
    anchor_.add_site(0, anchor_site_);
    mobile_.clear();
    mobile_.add_site(0, mobile_site_);
  }

  bool select(const Select& perturbed, System * system) override {
    Configuration * config = system->get_configuration();
    int particle_index = -1;
    if (perturbed.num_sites() > 0) {
      particle_index = perturbed.particle_indices().back();
    } else {
      // select random particle of correct type
      const int group_index = config->particle_type_to_group(particle_type());
      const int num = config->num_particles(group_index);
      if (num <= 0) return false;
      const int index = random_.uniform(0, num - 1);
      const SelectGroup& select = config->group_select(group_index);
      particle_index = select.particle_index(index);
    }
    mobile_.set_particle(0, particle_index);
    anchor_.set_particle(0, particle_index);
    mobile_.load_positions(config->particles());
    DEBUG("mobile: " << mobile_.str());
    DEBUG("anchor: " << anchor_.str());
    mobile_original_ = mobile_;
    return true;
  }

 private:
  int mobile_site_;
  int anchor_site_;
  Random random_;
};

inline std::shared_ptr<TrialSelectBond> MakeTrialSelectBond(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectBond>(args);
}

/**
  A random particle of given type is selected if previously perturbed sites
    are not available.
  Select a single angle from two given anchor sites and one mobile site.
  The mobile site is directly bonded to the first anchor site, and the second
    anchor site is bonded to the first mobile site.
  The angle is defined as: anchor2 - anchor1 - mobile.
 */
class TrialSelectAngle : public TrialSelectBond {
 public:
  TrialSelectAngle(
    /**
      anchor_site2 : index of second anchor site.
     */
    const argtype& args = argtype()) : TrialSelectBond(args) {
    Arguments args_(args);
    args_.dont_check();
    anchor_site2_ = args_.key("anchor_site2").integer();
  }

  // angle theta0 is added as a property
  // anchor is sized
  void precompute(System * system) override {
    TrialSelectBond::precompute(system);
    const Particle& part = system->configuration().particle_types().particle(particle_type());
    const int angle_type = part.angle(mobile_site(),
                                      anchor_site(),
                                      anchor_site2_).type();
    const Angle& angle = system->configuration().unique_types().particle(
      particle_type()).angle(angle_type);
    add_property("theta0", angle.property("theta0"));
    anchor_.add_site(0, anchor_site2_);
  }

 private:
  int anchor_site2_;
};

inline std::shared_ptr<TrialSelectAngle> MakeTrialSelectAngle(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_H_
