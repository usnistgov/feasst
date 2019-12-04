
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_H_

#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "system/include/select_list.h"

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
    Random * random) {
    ERROR("not implemented"); }

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(System * system);

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
    const System * system);

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

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<TrialSelect> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<TrialSelect> >& deserialize_map();
  std::shared_ptr<TrialSelect> deserialize(std::istream& istr);
  virtual ~TrialSelect() {}

 protected:
  std::string class_name_ = "TrialSelect";
  SelectList mobile_original_;
  SelectList mobile_;
  Select anchor_;

  void serialize_trial_select_(std::ostream& ostr) const;
  TrialSelect(std::istream& istr);

 private:
  int group_index_;
  int particle_type_;
  bool is_particle_type_set_ = false;
  bool is_ghost_;

  // not checkpointed
  double probability_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_H_
