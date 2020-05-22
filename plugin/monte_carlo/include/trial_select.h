
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_H_

#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "configuration/include/properties.h"
#include "configuration/include/select.h"
#include "system/include/system.h"

namespace feasst {

/**
  Select the mobile particles and sites that are to be perturbed via trials.
  Store the original position in mobile_original for reverting.
  Also store the 'anchor' particles and sites which are not mobile but may be
  required to complete the perturbation (e.g., bonds).
 */
class TrialSelect : public PropertiedEntity {
 public:
  /**
    args:
    - group_index: index of group definied within system (default: 0).
    - particle_type: type of particle in configuration (default: -1)
   */
  explicit TrialSelect(const argtype& args = argtype());

  /// Return the index of group for selection.
  int group_index() const { return group_index_; }

  /// Return the particle type.
  int particle_type() const;

  /// Perform upkeep before select.
  void before_select() { mobile_.reset_excluded_and_bond(); }

  /// Perform the selection as implemented in the derived class.
  /// Return false if the selection cannot be made. Otherwise, return true.
  virtual bool select(
    /// Perturbed is included to allow chaining of selection based on previous.
    const Select& perturbed,
    System * system,
    Random * random);

  /// Same as above but with an empty perturbed.
  bool sel(System * system, Random * random) {
    return select(empty_, system, random); }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(System * system);

  /// Return the mobile selection. These can change during the trial.
  const Select& mobile() const { return mobile_; }

  /// Return a pointer to the mobile selection.
  Select * get_mobile() { return &mobile_; }

  /// Set the mobile selection.
  void set_mobile(const Select& mobile) { mobile_ = mobile; }

  /// Return originally-seleted mobile. These do not change during trial.
  const Select& mobile_original() const { return mobile_original_; }

  /// Return the anchor selection.
  const Select& anchor() const { return anchor_; }

  /// Return anchor position.
  const Position& anchor_position(
    /// anchor index, not configuration index
    const int particle_index,
    /// anchor index
    const int site_index,
    const System& system);

  /// Set the state of the trial for the mobile select (e.g., old, move, add).
  /// See Select::trial_state
  void set_trial_state(const int state) { mobile_.set_trial_state(state); }

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
  void set_ghost(const bool ghost = true);

  /// Return true if selecting from ghost particles.
  bool is_ghost() const { return is_ghost_; }

  /// Return printable properties.
  const std::map<std::string, Accumulator>& printable() const { return printable_; }
  const Accumulator& printable(const std::string str) { return printable_[str]; }

  /// Return true if constraints are satisfied.
  virtual bool are_constraints_satisfied(const System& system) const {
    return true; }

  /// Return true if particle type is set.
  bool is_particle_type_set() const { return is_particle_type_set_; }

  /// Remove unphysical sites from mobile.
  void remove_unphysical_sites(const Configuration& config);

  /// Fast replace of a single particle in mobile.
  void replace_mobile(const Select& replacement, const int sel_part_index,
    const Configuration& config);

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<TrialSelect> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<TrialSelect> >& deserialize_map();
  std::shared_ptr<TrialSelect> deserialize(std::istream& istr);
  virtual ~TrialSelect() {}

 protected:
  std::string class_name_ = "TrialSelect";
  Arguments args_;
  Select mobile_original_;
  Select mobile_;
  Select anchor_;
  std::map<std::string, Accumulator> printable_;

  void serialize_trial_select_(std::ostream& ostr) const;
  TrialSelect(std::istream& istr);

  const EnergyMap& map_(const System& system,
    const NeighborCriteria& neighbor_criteria) const;

 private:
  int group_index_;
  int particle_type_;
  bool is_particle_type_set_ = false;
  bool is_ghost_;

  // not checkpointed
  double probability_;
  Select empty_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_H_
