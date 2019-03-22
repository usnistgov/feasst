
#ifndef FEASST_CORE_PERTURB_H_
#define FEASST_CORE_PERTURB_H_

#include "core/include/system.h"
#include "core/include/select_list.h"
#include "core/include/tunable.h"

namespace feasst {

/**
  Perturb the system (e.g., displace, add or delete particles).
  Importantly, these moves are reversible upon calling the revert function.
  HWH refactor: consider having perturb only work with configuration and not system.
 */
class Perturb {
 public:
  virtual void perturb(System * system) { ERROR("not implemented"); }

  /// Initialize some variables before each attempt.
  virtual void before_attempt() { revert_possible_ = false; }

  /// Revert the perturbation.
  virtual void revert() {
    ASSERT(!optimized_revert(), "nonoptimized revert requires system storage");
    *system_ = system_old_;
  }

  /// Return whether it is possible to revert.
  bool revert_possible() const { return revert_possible_; }

  /// Set whether it is possible to revert.
  void set_revert_possible(const bool revert = true) {
    revert_possible_ = revert;
  }

  /// Return the selection
  //virtual const Select& selection() const = 0;
  const SelectList& selection() const { return selection_; }

  /// Set the selection
  void set_selection(const SelectList& select) { selection_ = select; }

  void set_selection_state(const std::string state) {
    selection_.set_trial_state(state); }

  void select_random_particle(const int group_index, const Configuration& config) {
    selection_.random_particle(config, group_index);
    set_selection_state("old");
  }

  void select_last_particle_added(const Configuration * config) {
    selection_.last_particle_added(config); }

  void select_random_particle_of_type(const int type, Configuration * config,
    /// Don't load coordinates into selection by default
    const int load_coordinates = 0) {
    selection_.random_particle_of_type(type, config, load_coordinates);
  }

  /// Select all sites between two randomly selected sites in a randomly selected particle in group.
  void select_random_segment_in_particle(const int group_index, const Configuration& config, const int max_length = -1) {
    selection_.random_segment_in_particle(group_index, config, max_length);
    set_selection_state("old");
  }

  /// Select all sites between a random endpoint and a randomly selectioned site in a randomly selected particle in group.
  /// Note that the end point is always returned as the first site.
  /// Thus, when the end point is the last site, the site order is reversed.
  /// HWH note: copy of comment in SelectList
  void select_random_end_segment_in_particle(const int group_index, const Configuration& config, const int max_length = -1) {
    selection_.random_end_segment_in_particle(group_index, config, max_length);
    set_selection_state("old");
  }

  void set_tunable(const Tunable tunable) { tunable_ = tunable; }
  Tunable tunable() const { return tunable_; }
  void tune(const double value) { tunable_.tune(value); }
  void set_tune_min_and_max(const double min, const double max) {
    tunable_.set_min_and_max(min, max); }

  /// Set the anchors, which are not moved in the perturbation but are used
  /// for reference.
  void set_anchor(const SelectList& anchor) { anchor_ = anchor; }
  const SelectList& anchor() const { return anchor_; }

  virtual ~Perturb() {}

 protected:
  /* The following protected functions are only to be used by developers */

  // Before each perturbation, store the old system.
  void store_old(System * system) {
    system_ = system;
    if (!optimized_revert()) {
      system_old_ = *system;
    }
  }

  virtual bool optimized_revert() { return false; }

  System* system() { return system_; }

 private:
  System * system_;
  System system_old_;
  bool revert_possible_;
  SelectList selection_;
  SelectList anchor_;
  Tunable tunable_;
};

class PerturbOptRevert : public Perturb {
  bool optimized_revert() override { return true; }
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_H_
