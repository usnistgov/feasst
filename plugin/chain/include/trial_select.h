
#ifndef FEASST_CHAIN_TRIAL_SELECT_H_
#define FEASST_CHAIN_TRIAL_SELECT_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select a random segment.
// HWH optimzie by settting endpoints as anchors.
class TrialSelectSegment : public TrialSelectParticle {
 public:
  TrialSelectSegment(
    /**
      max_length : maximum length of selected segment. If -1 (default), then
        randomly select all possible lengths.
     */
    const argtype& args = argtype()) : TrialSelectParticle(args) {
    Arguments args_(args);
    args_.dont_check();
    max_length_ = args_.key("max_length").dflt("-1").integer();
  }

  int max_length() const { return max_length_; }

  /// Select all sites between two randomly selected sites in a randomly selected particle in group.
  void random_segment_in_particle(
      const Configuration& config,
      SelectPosition * select,
      Random * random,
      /// Set the maximum length of the segment.
      /// If -1 (default), consider all possible lengths.
      const int max_length = -1
    ) {
    random_particle(config, select, random);
    const int num_sites = select->num_sites();
    if (num_sites <= 1) {
      return; // HWH note this check prevents error/infinite loop below
    }

    // find two unequal sites
    int min = 0;
    int max = min;
    int attempt = 0;
    while (min == max) {
      min = random->uniform(0, num_sites - 1);
      if (max_length == -1) {
        max = random->uniform(0, num_sites - 1);
      } else {
        max = min + random->uniform(-max_length, max_length);
        if (max < 0) {
          max = 0;
        }
        if (max >= num_sites) {
          max = num_sites - 1;
        }
      }
      ++attempt;
      ASSERT(attempt < 1e3, "infinite loop");
    }

    // swap for meaningful min/max
    sort(&min, &max);

    // remove sites not in min/max, from highest to lowest
    select->remove_last_sites(num_sites - max - 1);
    select->remove_first_sites(min);
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    random_segment_in_particle(
      system->configuration(),
      &mobile_,
      random,
      max_length()
    );
    mobile_original_ = mobile_;
    return true;
  }

  virtual ~TrialSelectSegment() {}

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

  /// Select all sites between a random endpoint and a randomly selectioned site in a randomly selected particle in group.
  /// Return true if the endpoint is at the beginning.
  bool random_end_segment_in_particle(const Configuration& config,
      SelectPosition * select,
      Random * random,
      /// Set the maximum length of the segment.
      /// If -1 (default), consider all possible lengths.
      const int max_length = -1
      ) {
    random_particle(config, select, random);
    const int num_sites = select->num_sites();
    // HWH note this check prevents error/infinite loop below
    if (num_sites <= 1) {
      DEBUG("num sites(" << num_sites << ") not large enough");
      return false;
    }

    // select a random site
    int site = -1;
    bool is_endpoint_beginning;
    if (max_length == -1) {
      site = random->uniform(0, num_sites - 1);

      DEBUG("site " << site << " num " << num_sites);

      // randomly decide which endpoint to keep in selection
      if (site == 0) {
        is_endpoint_beginning = false;
      } else if (site == num_sites - 1) {
        is_endpoint_beginning = true;
      } else {
        if (random->coin_flip()) {
          is_endpoint_beginning = false;
        } else {
          is_endpoint_beginning = true;
        }
      }
    } else {
      ASSERT(max_length > 0, "max_length(" << max_length <<") should be >0 "
        << "or no segment will be selected");
      if (random->coin_flip()) {
        is_endpoint_beginning = false;
        site = random->uniform(num_sites - max_length, num_sites - 1);
      } else {
        is_endpoint_beginning = true;
        site = random->uniform(0, max_length - 1);
      }
    }

    DEBUG("beginning? " << is_endpoint_beginning);
    if (is_endpoint_beginning) {
      select->remove_last_sites(num_sites - site - 1);
    } else {
      select->remove_first_sites(site);
    }
    DEBUG("num " << num_sites << " indices " << select->str());
    return is_endpoint_beginning;
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    bool is_endpoint_beginning = random_end_segment_in_particle(
      system->configuration(),
      &mobile_,
      random,
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

  virtual ~TrialSelectEndSegment() {}
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

  virtual ~TrialSelectReptate() {}

 private:
  // temporary
  SelectList bonded_to_;
};

/// Select the perturbed site. Used for selection of first stage with growth
/// expected ensemble. Particular particle, not random particle.
class TrialSelectPerturbed : public TrialSelect {
 public:
  TrialSelectPerturbed(const argtype& args = argtype()) : TrialSelect(args) {}

  bool select(const Select& perturbed, System* system, Random * random) override {
    if (perturbed.num_sites() == 0) return false;
    mobile_.replace_particle(perturbed, 0, system->configuration());
    mobile_original_ = mobile_;
    return true;
  }

  virtual ~TrialSelectPerturbed() {}
};

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_H_
