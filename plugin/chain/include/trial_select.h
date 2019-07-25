
#ifndef FEASST_CHAIN_TRIAL_SELECT_H_
#define FEASST_CHAIN_TRIAL_SELECT_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select a random segment.
// HWH optimzie by settting endpoints as anchors.
class TrialSelectSegment : public TrialSelect {
 public:
  TrialSelectSegment(
    /**
      max_length : maximum length of selected segment. If -1 (default), then
        randomly select all possible lengths.
     */
    const argtype& args = argtype()) : TrialSelect(args) {
    Arguments args_(args);
    args_.dont_check();
    max_length_ = args_.key("max_length").dflt("-1").integer();
  }

  int max_length() const { return max_length_; }

  bool select(const Select& perturbed, System* system) override {
    mobile_.random_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length()
    );
    mobile_original_ = mobile_;
    return true;
  }

 private:
  int max_length_;
};

/// Select an end segment.
/// Set the anchor as the other end of the selection from the end point.
// HWH optimize, set anchor one site next from selection.
class TrialSelectEndSegment : public TrialSelectSegment {
 public:
  TrialSelectEndSegment(const argtype& args = argtype()) : TrialSelectSegment(args) {}
  void precompute(System * system) override {
    anchor_.clear();
    anchor_.add_site(0, 0);
  }

  bool select(const Select& perturbed, System* system) override {
    bool is_endpoint_beginning = mobile_.random_end_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length()
    );
    if (mobile_.num_sites() <= 0) {
      return false;
    }
    update_anchor(is_endpoint_beginning, system);
    mobile_original_ = mobile_;
    return true;
  }

  virtual void update_anchor(const bool is_endpoint_beginning,
    const System * system) {
    int select_index = -1;
    if (is_endpoint_beginning) {
      select_index = mobile_.num_sites() - 1;
    } else {
      select_index = 0;
    }
    DEBUG("is_endpoint_beginning, " << is_endpoint_beginning);
    DEBUG("site index " << select_index);
    anchor_.set_site(0, 0, mobile_.site_indices()[0][select_index]);
    anchor_.set_particle(0, mobile_.particle_indices()[0]);
  }
};

const argtype TrialSelectReptateArg_ = {{"max_length", "1"}};

/// Select a random end point for reptation.
class TrialSelectReptate : public TrialSelectEndSegment {
 public:
  TrialSelectReptate() : TrialSelectEndSegment(TrialSelectReptateArg_) {}

  void precompute(System * system) override {
    anchor_.clear();
    anchor_.add_site(0, 0);
    bonded_to_.clear();
    bonded_to_.add_site(0, 0);
  }

  void update_anchor(const bool is_endpoint_beginning,
    const System * system) override {
    const int particle_index = mobile_.particle_indices()[0];
    const Configuration& config = system->configuration();
    const Particle& particle = config.select_particle(particle_index);
    int anchor_index = 0;
    int site_bonded_to = particle.num_sites() - 2;
    DEBUG("is_endpoint_beginning " << is_endpoint_beginning);
    if (is_endpoint_beginning) {
      anchor_index = particle.num_sites() - 1;
      site_bonded_to = 1;
    }
    // for the old configuration, set the anchor to the old bond.
    anchor_.set_site(0, 0, anchor_index);
    anchor_.set_particle(0, particle_index);
    ASSERT(bonded_to_.replace_indices(particle_index, {site_bonded_to}),
      "bonded_to_ wasn't initialized to proper size on precompute");
  }

  void mid_stage() override {
    // exclude the anchor from interactions.
    // include interactions with site that use to be bonded
    mobile_.set_new_bond(anchor_);
    mobile_.set_old_bond(bonded_to_);
  }

 private:
  // temporary
  SelectList bonded_to_;
};

/// Select the perturbed site. Used for selection of first stage with growth
/// expected ensemble. Particular particle, not random particle.
class TrialSelectPerturbed : public TrialSelect {
 public:
  TrialSelectPerturbed(const argtype& args = argtype()) : TrialSelect(args) {}

  bool select(const Select& perturbed, System* system) override {
    if (perturbed.num_sites() == 0) return false;
    mobile_.replace_particle(perturbed, 0, system->configuration());
    mobile_original_ = mobile_;
    return true;
  }
};

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_H_
