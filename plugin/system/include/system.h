
#ifndef FEASST_SYSTEM_SYSTEM_H_
#define FEASST_SYSTEM_SYSTEM_H_

#include <vector>
#include <memory>
#include "configuration/include/configuration.h"
#include "system/include/potential_factory.h"

namespace feasst {

/**
  System is a facade design pattern in order to constrain and/or simplify
  the interface with multiple configurations and multiple lists of potentials.
  The typing and grouping of multiple configurations should be the same.
  There are three types of potentials.

  1. The first is without optimizations.

  2. The second is with optimizations, which may be periodically compared
     with the first.

  3. The remaining potentials are used for reference.
     For example, cheap energy calculations in dual-cut configurational bias.
     Also Mayer-sampling.
 */
class System {
 public:
  System() {}

  /** @name Configurations
    Store and retrieve a list of configurations. */
  //@{

  /// Add a configuration.
  void add(const Configuration& configuration);

  /// Return the configuration
  const Configuration& configuration(const int config = 0) const;

  /// Get the configuration as a pointer.
  /// This interface is to be avoided if possible.
  Configuration* get_configuration(const int config = 0);

  /// Return the dimensionality of the system.
  int dimension(const int config = 0) const;

  //@}
  /** @name Potentails
    Store and retrieve a list of potentials.
   */
  //@{

  /// Add a potential. By default, this is unoptimized.
  void add(const Potential& potential);

  /// Add an unoptimized potential.
  void add_to_unoptimized(const Potential& potential);

  /// Return the unoptimized potentials.
  const PotentialFactory& unoptimized() const { return unoptimized_; }

  /// Return an unoptimized potential.
  const Potential& potential(const int index) const {
    return unoptimized_.potentials()[index]; }

  /// Add an optimized potential.
  void add_to_optimized(const Potential& potential);

  /// Return the optimized potentials.
  const PotentialFactory& optimized() const { return optimized_; }

  /// Add a reference potential.
  void add_to_reference(const Potential& ref,
    /// Store different references by index.
    const int index = 0);

  /// Return a reference potential.
  const Potential& reference(const int ref, const int potential) const;

  /// Return the list of reference potentials.
  const std::vector<PotentialFactory> references() const { return references_; }

  //@}
  /** @name Energy
    Compute energies using a combination of configurations and potentials.
   */

  /// Precompute quantities for optimizations before calculation of energies.
  void precompute();

  /// Return the unoptimized energy. The following use optimized if available.
  double unoptimized_energy(const int config = 0);

  /// Return the energy of all.
  // HWH note:
  // for when bonded energies, etc come into play,
  // perhaps have energy of full system use reference potentials.
  // this way CB could distinguish external and internal interactions.
  double energy(const int config = 0);

  /// Return the energy of the selection.
  double energy(const Select& select, const int config = 0);

  /// Return the last computed energy.
  double stored_energy() const {
    return const_potentials_()->stored_energy(); }

  /// Return the profile of energies that were last computed.
  std::vector<double> stored_energy_profile() const {
    return const_potentials_()->stored_energy_profile(); }

  /// Return the reference energy.
  double reference_energy(const int ref = 0, const int config = 0) {
    return reference_(ref)->energy(&configurations_[0]); }

  /// Return the reference energy of the selection.
  double reference_energy(const Select& select, const int ref = 0,
    const int config = 0);

  /// Revert changes due to perturbations.
  void revert() { unoptimized_.revert(); }

  /// Set cache to load energy calculations.
  void load_cache(const bool load);

  /// Set cache to unload energy calclatuions.
  void unload_cache(const System& system);

  /// Serialize
  void serialize(std::ostream& sstr) const;

  /// Deserialize
  explicit System(std::istream& sstr);

 private:
  std::vector<Configuration> configurations_;
  PotentialFactory unoptimized_;
  PotentialFactory optimized_;
  bool is_optimized_ = false;
  std::vector<PotentialFactory> references_;

  PotentialFactory * reference_(const int index);
  PotentialFactory * potentials_();
  const PotentialFactory * const_potentials_() const;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_SYSTEM_H_
