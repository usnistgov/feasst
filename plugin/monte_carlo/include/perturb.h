
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
#include "math/include/accumulator.h"

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
  explicit Perturb(
    /**
      tunable_param : initial value of the tunable parameter (default: 0.1).
     */
    const argtype& args = argtype()) {
    Arguments args_(args);
    args_.dont_check();
    tunable_.set_value(args_.key("tunable_param").dflt("0.1").dble());
    before_select();
    set_probability();
  }

  /// Return the tunable parameter.
  const Tunable& tunable() const { return tunable_; }

  /// Set the minimum and maximum values of the tunable parameter.
  void set_tunable_min_and_max(const double min, const double max) {
    tunable_.set_min_and_max(min, max); }

  /// Set the value of the tunable parameter.
  void set_tunable(const double value) { tunable_.set_value(value); }

  /// Tune the parameter based on difference between target and actual.
  void tune(const double actual) { tunable_.tune(actual); }

  virtual void precompute(TrialSelect * select, System * system) {}

  /// Before perturbation, initialize some optimiation parameters.
  void before_select() {
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

  virtual std::string status_header() const {
    std::stringstream ss;
    if (tunable().is_enabled()) {
      ss << "tunable ";
    }
    return ss.str();
  }

  virtual std::string status() const {
    std::stringstream ss;
    if (tunable().is_enabled()) {
      ss << tunable_.value() << " ";
    }
    return ss.str();
  }

  virtual ~Perturb() {}

 protected:
  void disable_tunable_() { tunable_.disable(); }

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
      select->set_trial_state(0);
      return;
    }
    move(system, select);
    set_revert_possible(true, select);
    select->set_trial_state(1);
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

  virtual ~PerturbMove() {}
};

/**
  Translate the positions of the selection.
 */
class PerturbTranslate : public PerturbMove {
 public:
  PerturbTranslate(const argtype& args = argtype()) : PerturbMove(args) {}

  void precompute(TrialSelect * select, System * system) override {
    set_tunable_min_and_max(2*NEAR_ZERO,
      0.5*system->configuration().domain().max_side_length());
  }

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
    ASSERT(tunable().value() > NEAR_ZERO, "tunable(" << tunable().value()
      << ") is too small");
    move(trajectory_, system, select);
  }

  virtual ~PerturbTranslate() {}

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
  PerturbRotate(const argtype& args = argtype()) : PerturbMove(args) {
    set_tunable_min_and_max(2*NEAR_ZERO, 360.);
  }

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
    if (is_rotation_not_needed_(select, pivot)) return;
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

  virtual ~PerturbRotate() {}

 private:
  Random random_;

  const Position& piv_sel_(const Position& pivot, const TrialSelect * select) {
    if (pivot.dimension() == 0) {
      return select->mobile().particle_positions()[0];
    }
    return pivot;
  }

  // optimization
  // check to see if rotation is not necessary when pivot is equivalent to
  // the only rotated position.
  bool is_rotation_not_needed_(const TrialSelect * select,
      const Position& pivot) {
    const SelectList& rotated = select->mobile();
    if (rotated.num_particles() == 1) {
      if (static_cast<int>(rotated.site_indices()[0].size()) == 1) {
        if (rotated.site_positions()[0][0].is_equal(pivot)) {
          return true;
        }
      }
    }
    return false;
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
    DEBUG("anywhere: " << random_in_box_.str());
  }

  virtual ~PerturbAnywhere() {}

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
  PerturbAdd(const argtype& args = argtype()) : Perturb(args) {}

  //initialize ghost selection in TrialSelect?
  void precompute(TrialSelect * select, System * system) override {
    select->set_ghost(true);
  }

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
    DEBUG("is_position_held " << is_position_held);
    Configuration* config = system->get_configuration();
    config->revive(select->mobile());
    select->set_trial_state(2);

    // obtain probability
    const int particle_type = config->select_particle(
      select->mobile().particle_index(0)
    ).type();
    DEBUG("type " << particle_type);
    for (const Select& ghost : config->ghosts()) {
      DEBUG("ghost " << ghost.str());
    }
    select->set_probability(
      1./static_cast<double>(config->num_particles_of_type(particle_type)));

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
      DEBUG("nump " << system->configuration().num_particles());
      system->get_configuration()->remove_particles(revert_select()->mobile());
      system->revert();
    }
  }

  std::string status_header() const override {
    std::stringstream ss;
    return ss.str();
  }

  std::string status() const override {
    std::stringstream ss;
    return ss.str();
  }

  virtual ~PerturbAdd() {}

 private:
  PerturbAnywhere anywhere_;
  Select whole_particle_;

  // temporary
  Position empty_;
};

inline std::shared_ptr<PerturbAdd> MakePerturbAdd(const argtype& args = argtype()) {
  return std::make_shared<PerturbAdd>(args);
}

/**
  Remove a particle from the system.
 */
class PerturbRemove : public Perturb {
 public:
  PerturbRemove() {
    disable_tunable_();
  }

  void perturb(
      System * system,
      TrialSelect * select,
      const bool is_position_held = true
      ) override {
    select->set_trial_state(0);
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

  virtual ~PerturbRemove() {}

 private:
  PerturbAnywhere anywhere_;
};

/// Put first site in selection in a sphere about the first site in anchor.
/// HWH: could enable tuning and position w.r.t. previous bond placement
///      for higher acceptance probability
class PerturbDistance : public PerturbMove {
 public:
  PerturbDistance(const argtype& args = argtype())
    : PerturbMove(args) {
    disable_tunable_();
  }

  void precompute(TrialSelect * select, System * system) override {
    // determine the bond length
    // or input the bond length
    if (select->has_property("bond_length")) {
      distance_ = select->property("bond_length");
    } else {
      WARN("using default distance (typically for reptation): " << distance_);
    }
  }

  double distance() const { return distance_; }

  void move(System * system,
      TrialSelect * select) override {
    SelectList * mobile = select->get_mobile();
    Position * site = mobile->get_site_position(0, 0);
    DEBUG("mobile " << mobile->str());
    DEBUG("old pos " << site->str());
    random_.unit_sphere_surface(site);
    site->multiply(distance_);
    site->add(select->anchor_position(0, 0, system));
    DEBUG("new pos " << site->str());
    system->get_configuration()->update_positions(select->mobile());
  }

  virtual ~PerturbDistance() {}

 private:
  double distance_ = 1.;

  // temporary
  Random random_;
};

/// Put first site in selection, i, in a sphere about the first site in anchor,
///  j, and at an angle i,j,k (vertex: j) about the second site in anchor, j.
class PerturbDistanceAndAngle : public PerturbDistance {
 public:
  PerturbDistanceAndAngle(const argtype& args = argtype())
    : PerturbDistance(args) {}

  void precompute(TrialSelect * select, System * system) override {
    PerturbDistance::precompute(select, system);
    angle_ = select->property("theta0");
    origin_.set_to_origin(system->configuration().dimension());
    rjk_ = origin_;
  }

  void move(System * system,
      TrialSelect * select) override {
    SelectList * mobile = select->get_mobile();
    Position * site = mobile->get_site_position(0, 0);
    DEBUG("mobile " << mobile->str());
    DEBUG("old pos " << site->str());

    // set site to the vector |r_j - r_k| and store this unit vector
    const Position& rj = select->anchor_position(0, 0, system);
    const Position& rk = select->anchor_position(0, 1, system);
    *site = rj;
    site->subtract(rk);
    rjk_ = *site;
    DEBUG("rjk " << rjk_.str());
    //rjk_.normalize();

    // rotate site by (PI-theta0) about vector orthogonal to r_jk
    orthogonal_jk_.orthogonal(*site);
    DEBUG("ortho " << orthogonal_jk_.str());
    rot_mat_.axis_angle(orthogonal_jk_, PI - angle_);
    DEBUG("site == rj: " << site->str());
    rot_mat_.rotate(origin_, site);
    DEBUG("site rotated to angle: " << site->str());

    // randomly spin site about rjk.
    rot_mat_.axis_angle(rjk_, 2*PI*random_.uniform());
    rot_mat_.rotate(origin_, site);

    site->add(rj);  // return frame of reference

    DEBUG("new pos " << site->str());
    system->get_configuration()->update_positions(select->mobile());
  }

  virtual ~PerturbDistanceAndAngle() {}

 private:
  double angle_;

  // temporary
  Random random_;
  Position rjk_;
  Position orthogonal_jk_;
  Position origin_;
  RotationMatrix rot_mat_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_H_
