
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_H_

#include <memory>
#include <map>
#include <string>

namespace feasst {

class Accumulator;
class Configuration;
class EnergyMap;
class Position;
class Properties;
class Random;
class Select;
class System;

typedef std::map<std::string, std::string> argtype;

/**
  Select the mobile particles and sites that are to be perturbed via trials.
  Store the original position in mobile_original for reverting.
  Also store the 'anchor' particles and sites which are not mobile but may be
  required to complete the perturbation (e.g., bonds).
 */
class TrialSelect {
 public:
  //@{
  /** @name Arguments
    - group_index: index of group defined within system (default: 0).
    - group: name of group defined within system (default: "").
    - particle_type: particle type name in Configuration (default: not set)
    - config: name of Configuration (default: 0)
   */
  explicit TrialSelect(argtype args = argtype());
  explicit TrialSelect(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the index of group for selection.
  int group_index() const { return group_index_; }

  /// Return the particle type.
  int particle_type() const;

  /// Return the name of the particle type.
  const std::string& particle_type_name() const { return particle_type_name_; }

  /// Return the index of the configuration.
  int configuration_index() const { return configuration_index_; }

  /// Given the system, return the configuration.
  const Configuration& configuration(const System& system) const;

  /// Given the system pointer, return the configuration pointer.
  Configuration * get_configuration(System * system) const;

  /// Set the configuration index.
  void set_configuration_index(const int config);

  /// Set the config name.
  void set_config(const std::string& config);

  /// Return the config name.
  const std::string& config() { return config_; }

  /// Perform upkeep before select.
  void before_select();

  /// Perform the selection as implemented in the derived class.
  /// Return false if the selection cannot be made. Otherwise, return true.
  virtual bool select(
    /// Perturbed is included to allow chaining of selection based on previous.
    const Select& perturbed,
    System * system,
    Random * random,
    TrialSelect * previous_select);

  /// Same as above but with an empty perturbed.
  bool sel(System * system, Random * random);

  /// Precompute quantities before simulation for optimization.
  virtual void precompute(System * system);

  /// Return the mobile selection. These can change during the trial.
  const Select& mobile() const;

  // Return a pointer to the mobile selection.
  Select * get_mobile();

  /// Set the mobile selection.
  void set_mobile(const Select& mobile);

  /// Return originally-seleted mobile. These do not change during trial.
  const Select& mobile_original() const;

  /// Set the original mobile, including Euler if anisotropic.
  void set_mobile_original(const System * system);

  /// Return the anchor selection.
  const Select& anchor() const;

  // Return a pointer to the anchor selection.
  Select * get_anchor();

  /// Return anchor position.
  const Position& anchor_position(
    /// anchor index, not configuration index
    const int particle_index,
    /// anchor index
    const int site_index,
    const System& system) const;

  /// Set the state of the trial for the mobile select (e.g., old, move, add).
  /// See Select::trial_state
  void set_trial_state(const int state);

  /// Reset the mobile selection to the original.
  void reset_mobile();

  /// Return the probability of the selection. For example, if a random particle
  /// type is selected, then the probability is the inverse of the number of
  /// particles of that type.
  double probability() const { return probability_; }

  /// Call between rosenbluth calculation of old and new configuration.
  virtual void mid_stage() {}

  /// Select from ghost particles.
  void set_ghost(const bool ghost = true);

  /// Return true if selecting from ghost particles.
  bool is_ghost() const { return is_ghost_; }

  /// Return printable properties.
  const std::map<std::string, std::shared_ptr<Accumulator> >& printable() const;
  const Accumulator& printable(const std::string str) const;

  /// Return true if constraints are satisfied.
  virtual bool are_constraints_satisfied(const int old,
    const System& system) const;

  /// Return true if particle type is set.
  bool is_particle_type_set() const { return is_particle_type_set_; }

  /// Remove unphysical sites from mobile.
  void remove_unphysical_sites(const Configuration& config);

  /// Fast replace of a single particle in mobile.
  void replace_mobile(const Select& replacement, const int sel_part_index,
    const Configuration& config);

  // HWH used in ComputeAddAVBDivalent
  const EnergyMap& map_(const System& system, const int neighbor_index) const;

  /// Return the property value by name.
  double property(const std::string name) const;

  /// Return true if entity has property of name.
  bool has_property(const std::string name) const;

  /// Add a property, or set its value if name already exists.
  void add_or_set_property(const std::string name, const double value);

  /// Return true if all sites in mobile are isotropic.
  bool is_isotropic(const System * system) const;

  // HWH hackish interface to exclude bond energies from lnpmet
  void zero_exclude_energy() { exclude_energy_ = 0.; }
  void add_exclude_energy(const double energy);
  double exclude_energy() const { return exclude_energy_; }

  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<TrialSelect> create(std::istream& istr) const;
  std::map<std::string, std::shared_ptr<TrialSelect> >& deserialize_map();
  std::shared_ptr<TrialSelect> deserialize(std::istream& istr);
  virtual ~TrialSelect() {}

  //@}
 protected:
  std::string class_name_ = "TrialSelect";
  std::shared_ptr<Select> mobile_original_;
  std::shared_ptr<Select> mobile_;
  std::shared_ptr<Select> anchor_;
  std::map<std::string, std::shared_ptr<Accumulator> > printable_;

  /// Set the probability of selection.
  void set_probability_(const double prob = 1) { probability_ = prob; }

  void serialize_trial_select_(std::ostream& ostr) const;
  explicit TrialSelect(std::istream& istr);

 private:
  int group_index_;
  std::string group_;
  int particle_type_;
  std::string particle_type_name_;
  int configuration_index_;
  std::string config_;
  bool is_particle_type_set_ = false;
  bool is_ghost_;
  std::shared_ptr<Properties> properties_;
  int aniso_index_ = -1;

  // not checkpointed
  double probability_;
  std::shared_ptr<Select> empty_;
  double exclude_energy_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_H_
