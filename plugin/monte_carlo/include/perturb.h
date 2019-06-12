
#ifndef FEASST_MONTE_CARLO_PERTURB_H_
#define FEASST_MONTE_CARLO_PERTURB_H_

#include <vector>
#include <numeric>
#include <string>
#include <memory>
#include "system/include/system.h"
#include "monte_carlo/include/tunable.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/**
  Perturbations to the system include the following types:

    1. move (single particle, regrowth, cluster moves).
    2. add/remove (grand canonical).
    3. change type (swap, growth expanded).
    4. change parameter (temperature, criteria swap bound).

  Perturbations may be followed up by one of the two following options:

    1. Revert the perturbation. For example, if a particle is moved, then
       restore the positions of that particle to the original.

    2. Finalize the perturbation. For example, if a selected particle is
       accepted for deletion, the particle is deleted in the finalize step.
 */
class Perturb {
 public:
  Perturb(
    /**
      tunable_param : initial value of the tunable parameter (default: 0.1).
     */
    const argtype& args = argtype()) {
    Arguments args_(args);
    args_.dont_check();
    tunable_.set_value(args_.key("tunable_param").dflt("0.1").dble());
    before_perturb();
    set_probability();
  }

  /// Return the tunable parameter.
  Tunable tunable() const { return tunable_; }

  /// Set the value of the tunable parameter.
  void set_tunable(const double value) { tunable_.set_value(value); }

  /// Tune the parameter based on difference between target and actual.
  void tune(const double actual) { tunable_.tune(actual); }

  virtual void precompute(const TrialSelect * select, System * system) {}

  /// Before perturbation, initialize some optimiation parameters.
  void before_perturb() {
    revert_possible_ = false;
    finalize_possible_  = false;
  }

  /// Perturb the selection of the system.
  virtual void perturb(
    System * system,
    TrialSelect * select,
    /// If position is held, all but perform the actual perturbation, as typically
    /// required for calculations of old configurations and Rosenbluth factors.
    const bool is_position_held = false
    ) { ERROR("not implemented"); }

  /// Return whether it is possible to revert.
  bool revert_possible() const { return revert_possible_; }

  /// Set whether it is possible to revert.
  void set_revert_possible(const bool revert_possible,
    /// If possible, store the selection.
    TrialSelect * revert_select) {
    revert_possible_ = revert_possible;
    revert_select_ = revert_select;
  }

  /// Return the stored revert selection.
  const TrialSelect * revert_select() const { return revert_select_; }

  /// Revert the system for an unaccepted perturbation.
  virtual void revert(System * system) {
    if (revert_possible_) {
      ERROR("not implemented");
    }
  }

  /// Return whether it is possible to finalize.
  bool finalize_possible() const { return finalize_possible_; }

  /// Set whether it is possible to finalize.
  void set_finalize_possible(const bool finalize_possible,
    /// If possible, store the selection.
    TrialSelect * finalize_select) {
    finalize_possible_ = finalize_possible;
    finalize_select_ = finalize_select;
    DEBUG("finalize possible: " << finalize_possible_);
  }

  /// Return the stored finalize selection.
  const TrialSelect * finalize_select() const { return finalize_select_; }

  /// Finalize the system for an accepted perturbation.
  virtual void finalize(System * system) {
    if (finalize_possible_) {
      ERROR("not implemented");
    }
  }

  /// Return the probability.
  double probability() const { return probability_; }

  /// Set the probability.
  void set_probability(const double prob = 1) { probability_ = prob; }

 private:
  Tunable tunable_;

  // optimzation or temporary object
  bool revert_possible_, finalize_possible_;
  TrialSelect * revert_select_;
  TrialSelect * finalize_select_;
  double probability_;
};

/**
  Only perturb the positions of the particles and/or sites.
 */
class PerturbMove : public Perturb {
 public:
  PerturbMove(const argtype& args = argtype()) : Perturb(args) {}

  /// Move the selection of the system.
  virtual void move(System * system, TrialSelect * select) = 0;

  // The perturbation move is simplified such that the move of the selection of
  // the system is all that remains to be implemented.
  void perturb(
      System * system,
      TrialSelect * select,
      const bool is_position_held = false
      ) override {
    if (is_position_held) {
      select->set_trial_state("old");
      return;
    }
    move(system, select);
    set_revert_possible(true, select);
    select->set_trial_state("move");
  }

  /// For perturbations that only move particles and/or sites, the revert step
  /// is the same for all. Simply put the original positions back.
  void revert(System * system) override {
    if (revert_possible()) {
      Configuration* config = system->get_configuration();
      config->update_positions(revert_select()->mobile_original(),
        // don't wrap if reverting
        false);
      system->revert();
    }
  }
};

/**
  Translate the positions of the selection.
 */
class PerturbTranslate : public PerturbMove {
 public:
  PerturbTranslate(const argtype& args = argtype()) : PerturbMove(args) {}

  /// Change the position in the selection given a trajectory.
  void update_selection(const Position& trajectory,
      TrialSelect * select) {
    SelectList * displaced = select->get_mobile();
    for (int select_index = 0;
         select_index < displaced->num_particles();
         ++select_index) {
      Position displaced_part(displaced->particle_positions()[select_index]);
      displaced_part.add(trajectory);
      displaced->set_particle_position(select_index, displaced_part);
      for (int site = 0;
           site < static_cast<int>(displaced->site_indices(select_index).size());
           ++site) {
        Position displaced_site(displaced->site_positions()[select_index][site]);
        displaced_site.add(trajectory);
        displaced->set_site_position(select_index, site, displaced_site);
      }
    }
  }

  /// Move the selected particles given a trajectory.
  void move(
      const Position& trajectory,
      System * system,
      TrialSelect * select) {
    update_selection(trajectory, select);
    system->get_configuration()->update_positions(select->mobile());
  }

  /// Move the selected particles using the tuning parameter.
  void move(System * system, TrialSelect * select) override {
    random_.position_in_cube(
      system->dimension(),
      tunable().value(),
      &trajectory_
    );
    DEBUG("max move " << tunable().value());
    ASSERT(tunable().value() > NEAR_ZERO, "tunable is too small");
    move(trajectory_, system, select);
  }

 private:
  // optimization or temporary objects
  Random random_;
  Position trajectory_;
};

/**
  Rotate the positions of the selection.
 */
class PerturbRotate : public PerturbMove {
 public:
  PerturbRotate(const argtype& args = argtype()) : PerturbMove(args) {}

  /// Change the position in the selection given a pivot and rotation matrix.
  void update_selection(const Position& pivot,
      const RotationMatrix& rotation,
      TrialSelect * select,
      /// Rotate particle positions (default). Otherwise, do not.
      const bool rotate_particle_position = true) {
    SelectList * rotated = select->get_mobile();
    for (int select_index = 0;
         select_index < rotated->num_particles();
         ++select_index) {
      // rotate site positions
      for (int site = 0;
           site < static_cast<int>(rotated->site_indices(select_index).size());
           ++site) {
        Position position = rotated->site_positions()[select_index][site];
        rotation.rotate(pivot, &position);
        rotated->set_site_position(select_index, site, position);
      }

      // rotate or recenter particle positions
      if (rotate_particle_position) {
        Position position = rotated->particle_positions()[select_index];
        rotation.rotate(pivot, &position);
        rotated->set_particle_position(select_index, position);
      }
    }
  }

  /// Change the position of the selection given a pivot and rotation matrix.
  void move(const Position& pivot,
      const RotationMatrix& rotation,
      System * system,
      TrialSelect * select,
      /// Rotate particle positions (default). Otherwise, do not.
      const bool rotate_particle_position = true) {
    update_selection(pivot, rotation, select, rotate_particle_position);
    system->get_configuration()->update_positions(select->mobile());
  }

  /// Rotate the selected particles using the tuning parameter.
  /// Set the pivot using the first particle position, and also
  /// rotate the particle positions.
  void move(System * system,
      TrialSelect * select) override {
    ASSERT(select->mobile().num_sites() > 0, "selection error");
    const Position& pivot = select->mobile().particle_positions()[0];
    move(system, select, pivot, true);
  }

  /// Rotate the selected particles using the tuning parameter.
  void move(System * system,
      TrialSelect * select,
      /// If pivot is empty, use first particle position.
      const Position& pivot,
      /// Rotate particle positions if true. Otherwise, do not.
      const bool rotate_particle_position) {
    const double max_angle = tunable().value();
    ASSERT(std::abs(max_angle) > NEAR_ZERO, "max angle is too small");
    const Position& piv_sel = piv_sel_(pivot, select);
    move(piv_sel,
      random_.rotation(piv_sel.dimension(), max_angle),
      system,
      select,
      rotate_particle_position
    );
  }

  Random * random() { return &random_; }

 private:
  Random random_;

  const Position& piv_sel_(const Position& pivot, const TrialSelect * select) {
    if (pivot.dimension() == 0) {
      return select->mobile().particle_positions()[0];
    }
    return pivot;
  }
};

/// Rigidly move rigidly anywhere in the box with any orientation.
class PerturbAnywhere : public PerturbMove {
 public:
  PerturbAnywhere() {
    rotate_.set_tunable(180.);
  }

  void set_position(const Position& center, System * system, TrialSelect * select) {
    translate_.move(center, system, select);
  }

  void move(System * system, TrialSelect * select) override {
    ASSERT(std::abs(rotate_.tunable().value() - 180.) < NEAR_ZERO,
      "rotation tunable should be 180");
    rotate_.move(system, select);
    system->configuration().domain().random_position(&random_in_box_, &random_);
    set_position(random_in_box_, system, select);
  }

 private:
  PerturbTranslate translate_;
  PerturbRotate rotate_;

  // optimization or temporary
  Position random_in_box_;
  Random random_;
};

/**
  Add a particle to the system.
 */
// HWH optimize -> update cell list in finalize?
class PerturbAdd : public Perturb {
 public:
  PerturbAdd(
    /**
      particle_type: type of particle in configuration to add (default: 0).
     */
    const argtype& args = argtype()) : Perturb(args) {
    Arguments args_(args);
    args_.dont_check();
    particle_type_ = args_.key("particle_type").dflt("0").integer();
  }

  /// Return the particle type to add.
  int particle_type() const { return particle_type_; }

  void perturb(
      System * system,
      TrialSelect * select,
      const bool is_position_held = false
      ) override {
    add(system, select, empty_, is_position_held);
  }

  void add(
    System * system,
    TrialSelect * select,
    /// place particle anywhere if center is of zero dimension.
    const Position& center,
    const bool is_position_held = false
  ) {
    Configuration* config = system->get_configuration();
    config->add_particle_of_type(particle_type_);
    SelectList * added = select->get_mobile();
    added->last_particle_added(config);
    added->set_trial_state("add");
    if (center.dimension() == 0) {
      anywhere_.perturb(system, select, is_position_held);
    } else {
      anywhere_.set_position(center, system, select);
    }
    set_revert_possible(true, select);
  }

  void revert(System * system) override {
    if (revert_possible()) {
      DEBUG(revert_select()->mobile().str());
      system->get_configuration()->remove_particles(revert_select()->mobile());
      system->revert();
    }
  }

 private:
  PerturbAnywhere anywhere_;
  int particle_type_;

  // temporary
  Position empty_;
};

/**
  Remove a particle from the system.
 */
class PerturbRemove : public Perturb {
 public:
  PerturbRemove() {}

  void perturb(
      System * system,
      TrialSelect * select,
      const bool is_position_held = true
      ) override {
    select->set_trial_state("old");
    set_finalize_possible(true, select);

    if (is_position_held) {
      anywhere_.set_revert_possible(false, NULL);
    } else {
      anywhere_.perturb(system, select, is_position_held);
      set_revert_possible(true, select);
    }
  }

  void finalize(System * system) override {
    if (finalize_possible()) {
      system->get_configuration()->remove_particles(finalize_select()->mobile());
      // system->revert();
    }
  }

  void revert(System * system) override {
    if (revert_possible()) {
      anywhere_.revert(system);
    }
  }

 private:
  PerturbAnywhere anywhere_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_H_
