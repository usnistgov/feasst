
#ifndef FEASST_SYSTEM_SYSTEM_H_
#define FEASST_SYSTEM_SYSTEM_H_

#include <string>
#include <map>
#include <vector>
#include <memory>
#include "system/include/potential.h"
#include "system/include/potential_factory.h"

namespace feasst {

class BondVisitor;
class Configuration;
class NeighborCriteria;
class ThermoParams;

typedef std::map<std::string, std::string> argtype;

/**
  System is a facade design pattern in order to constrain and/or simplify
  the interface with multiple configurations and multiple lists of potentials.
  The typing and grouping of multiple configurations should be the same, but
  the domain and particle which physically exist may be different.
  There are three types of potentials.

  1. The first is without optimizations.

  2. The second is with optimizations, which may be periodically compared
     with the first.

  3. The remaining potentials are used for reference.
     For example, cheap energy calculations in dual-cut configurational bias.
     Also Mayer-sampling.

  Finally, the System also contains the ThermoParams.
 */
class System {
 public:
  /// Empty constructor
  System() {}

  // HWH make configuration a shared pointer so can be extended?
  // Construct with a Configuration.
  // System(std::make_shared<Configuration> config) { add(*config); }

  /** @name Configurations
    Store and retrieve a list of configurations. */
  //@{

  /// Add a configuration.
  void add(std::shared_ptr<Configuration> configuration);

  /// Return the number of configurations.
  int num_configurations() const;

  /// Return the configuration
  const Configuration& configuration(const int config = 0) const;

  // Get the configuration as a pointer.
  // This interface is to be avoided if possible.
  Configuration* get_configuration(const int config = 0);

  /// Return the index of the Configuration based on the name.
  int configuration_index(const std::string& name) const;

  /// Return the configuration based on the name.
  const Configuration& configuration(const std::string& name) const;

  /// Return the configuration pointer based on the name.
  Configuration * configuration(const std::string& name);

  /// Return the dimensionality of the system.
  int dimension(const int config = 0) const;

  //@}
  /** @name Potentails
    Store and retrieve a list of potentials.
   */
  //@{

  /// Add a potential. By default, the potential is considered unoptimized.
  void add(std::shared_ptr<Potential> potential);

  /// Set an unoptimized potential.
  void set_unoptimized(const int index, std::shared_ptr<Potential> potential,
    const int config = 0);

  /// Add an unoptimized potential.
  void add_to_unoptimized(std::shared_ptr<Potential> potential,
    const int config = 0);

  /// Return the unoptimized potentials.
  const PotentialFactory& unoptimized(const int config = 0) const {
    return unoptimized_[config]; }

  /// Return an unoptimized potential.
  const Potential& potential(const int index, const int config = 0) const {
    return unoptimized_[config].potential(index); }

  // Return an unoptimized potential.
  Potential * get_potential(const int index, const int config = 0) {
    return unoptimized_[config].get_potential(index); }

  /// Add an optimized potential.
  void add_to_optimized(std::shared_ptr<Potential> potential,
    const int config = 0);

  /// Return the optimized potentials.
  const PotentialFactory& optimized(const int config = 0) const {
    return optimized_[config]; }

  /// Add a reference potential.
  void add_to_reference(std::shared_ptr<Potential> ref,
    /// Store different references by index.
    const int index = 0,
    std::string name = "");

  /// Return the number of reference potentials.
  int num_references(const int config = 0) const;

  /// Return a reference potential.
  const Potential& reference(const int ref, const int potential,
    const int config = 0) const;

  /// Return the list of reference potentials.
  const std::vector<std::vector<PotentialFactory> > references() const {
    return references_; }

  /// Return the index of the reference potential based on the name and the
  /// Configuration.
  /// If the name is not found, return the number of current reference
  /// potentials.
  int reference_index(const int config, const std::string& name) const;

  /// Return a constant reference to the full potentials.
  const PotentialFactory& potentials(const int config = 0) const;

  /// Remove optimization when overlap is detected, which is default.
  void remove_opt_overlap();

  //@}
  /** @name Neighbor Criteria
    Define and store various criteria used for defining neighbors
   */
  //@{

  /// Add a NeighborCriteria.
  void add(std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const int config = 0);

  /// Return a NeighborCriteria by index in order added.
  const NeighborCriteria& neighbor_criteria(const int index,
    const int config) const;

  /// Return a NeighborCriteria by index in order added.
  const std::vector<std::shared_ptr<NeighborCriteria> >& neighbor_criteria(
    const int config) const;

  // Return a NeighborCriteria by index in order added.
  NeighborCriteria * get_neighbor_criteria(const int index,
                                           const int config);

  //@}
  /** @name Energy
    Compute energies using a combination of configurations and potentials.
   */
  //@{

  /// Precompute quantities for optimizations before calculation of energies.
  void precompute();

  /// Return the unoptimized energy. The following use optimized if available.
  double unoptimized_energy(const int config);

  /// Return the energy of all.
  double energy(const int config = 0);
  // HWH note:
  // for when bonded energies, etc come into play,
  // perhaps have energy of full system use reference potentials.
  // this way CB could distinguish external and internal interactions.

  /// Return the energy of the selection.
  /// But do not finalize this energy (e.g., Ewald, neighbors, etc).
  double perturbed_energy(const Select& select, const int config = 0);

  /// Return the last computed energy.
  double stored_energy(const int config = 0) const {
    return potentials(config).stored_energy(); }

  /// Return the profile of energies that were last computed.
  std::vector<double> stored_energy_profile(const int config = 0) const {
    return potentials(config).stored_energy_profile(); }

  /// Return the reference energy.
  double reference_energy(const int ref = 0, const int config = 0);

  /// Return the reference energy of the selection.
  double reference_energy(const Select& select,
    const int ref = 0,
    const int config = 0);

  /// Initialize and return total energy.
  double initialize(const int config = 0);

  //@}
  /** @name ThermoParams
    Store and retreive the thermodynamic parameters such as temperature, pressure, etc.
   */
  //@{

  /// Set the thermodynamic parameters.
  void set(std::shared_ptr<ThermoParams> thermo_params);

  /// Return the thermodynamic parameters.
  const ThermoParams& thermo_params() const;

  // Same as above, but as a constant pointer.
  const ThermoParams * thermo_params_ptr_() const;

  // Same as above, but as a pointer.
  ThermoParams * get_thermo_params() { return thermo_params_.get(); }

  /// Set the inverse temperature, \f$\beta\f$.
  void set_beta(const double beta);

  //@}
  // Other functions:

  /**
    Change the volume.

    args:
    - configuration: index of configuration (default: 0)
    - see Configuration for remaining arguments.
   */
  void change_volume(const double delta_volume, argtype args = argtype());
  void change_volume(const double delta_volume, argtype * args);

  /// Return the previous delta_volume.
  double delta_volume_previous() const { return delta_volume_previous_; }
//  /**
//    Change the volume in the opposite amount as the last volume change.
//    This is implemented to keep the total volume constant for Gibbs ensemble.
//    Return the actual volume change.
//
//    args:
//    - Same as change_volume above.
//   */
//  double constrained_volume_change(argtype * args);
//  double constrained_volume_change(argtype args = argtype());

  /// Return the total volume of all Configurations.
  double total_volume() const;

  /// Revert changes due to energy computation of perturbations.
  void revert(const Select& select, const int config = 0);

  /// Finalize changes due to energy computation of perturbations.
  void finalize(const Select& select, const int config = 0);
  void finalize(const int config = 0);

  /// Set cache to load energy calculations.
  void load_cache(const bool load);

  /// Set cache to unload energy calclatuions.
  void unload_cache(const System& system);

  void synchronize_(const System& system, const std::vector<std::shared_ptr<Select> >& perturbed);

  /// Return the header of the status for periodic output.
  std::string status_header() const;

  /// Return the brief status for periodic output.
  std::string status() const;

  /// Run checks.
  void check(const int config = 0) const;

  /// Serialize
  void serialize(std::ostream& sstr) const;

  /// Deserialize
  explicit System(std::istream& sstr);

 private:
  // The first index is the config, and the second (if available) is ref
  std::vector<std::shared_ptr<Configuration> > configurations_;
  std::vector<std::shared_ptr<BondVisitor> > bonds_;
  std::vector<PotentialFactory> unoptimized_;
  std::vector<PotentialFactory> optimized_;
  bool is_optimized_ = false;
  std::vector<std::vector<PotentialFactory> > references_;
  std::shared_ptr<ThermoParams> thermo_params_;

  // temporary variable, not needed for serialization
  // In order to finalize or restart the correct reference potential utilized
  // in a trial, this temporarily stores that reference potential index.
  int ref_used_last_ = -1;
  double delta_volume_previous_ = 1e30;  // implemented for Gibbs ensemble.

  PotentialFactory * reference_(const int index, const int config);
  PotentialFactory * potentials_(const int config);

  /// Return the index of the reference potential based on the name and the
  /// Configuration.
  /// If the name is not found, return the number of current reference
  /// potentials.
  int reference_index_(const int config, const std::string& name) const;
};

inline std::shared_ptr<System> MakeSystem() {
  return std::make_shared<System>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_SYSTEM_H_
