
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

  void select(System* system) override {
    mobile_.random_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length()
    );
    mobile_original_ = mobile_;
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

  void select(System* system) override {
    bool is_endpoint_beginning = mobile_.random_end_segment_in_particle(
      group_index(),
      system->configuration(),
      max_length()
    );
    if (mobile_.num_sites() <= 0) {
      return;
    }
    mobile_original_ = mobile_;
    update_anchor(is_endpoint_beginning, system);
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
    anchor_.set_site_position(0, 0, mobile_.site_positions()[0][select_index]);
    DEBUG("anchor pos " << anchor_.site_positions()[0][0].str());
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
    int site_index = 0;
    int site_bonded_to = 1;
    if (is_endpoint_beginning) {
      site_index = particle.num_sites() - 1;
      site_bonded_to = site_index - 1;
    }
    anchor_.set_site_position(0, 0, particle.site(site_index).position());

    // exclude the anchor from interactions.
    mobile_.set_new_bond(anchor_);

    // include interactions with site that use to be bonded
    ASSERT(bonded_to_.replace_indices(particle_index, {site_bonded_to}),
      "bonded_to_ wasn't initialized to proper size on precompute");
    mobile_.set_old_bond(bonded_to_);
  }

 private:
  // temporary
  SelectList bonded_to_;
};

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_H_
