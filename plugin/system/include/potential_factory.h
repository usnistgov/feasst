
#ifndef FEASST_SYSTEM_POTENTIAL_FACTORY_H_
#define FEASST_SYSTEM_POTENTIAL_FACTORY_H_

#include <vector>
#include <string>
#include <sstream>
#include "system/include/potential.h"
// #include "utils/include/timer.h"

namespace feasst {

/**
  A list of potentials.
 */
class PotentialFactory {
 public:
  PotentialFactory() {}

  /// Add potential.
  void add(const Potential& potential);

  /// Set a potential.
  void set(const int index, const Potential& potential);

  /// Return the potentials.
  const std::vector<Potential>& potentials() const { return potentials_; }

  /// Return the number of potentials.
  int num() const { return static_cast<int>(potentials_.size()); }

  /// Precompute quantities for optimizations before calculation of energies.
  void precompute(Configuration * config);

  /// Precompute a particular potential by index.
  void precompute(const int index, Configuration * config);

  /// Compute the energy of the given configuration.
  double energy(Configuration * config);

  /// Compute the energy of the selection in the configuration.
  double energy(const Select& select, Configuration * config);

  /// Return the profile of energies that were last computed.
  std::vector<double> stored_energy_profile() const;

  /// Return the last computed value of the energy.
  double stored_energy() const;

  /// Return a human-readable status.
  std::string str() const;

  /// Revert any changes to the configuration due to the last energy computation
  void revert();

  /// Finalize changes due to perturbations.
  void finalize();

  /// Set cache to load energy calculations.
  void load_cache(const bool load);

  /// Set cache to unload energy calclatuions.
  void unload_cache(const PotentialFactory& factory);

  /// Serialize.
  void serialize(std::ostream& sstr) const;

  /// Deserialize.
  explicit PotentialFactory(std::istream& sstr);

 private:
  std::vector<Potential> potentials_;
//  Timer timer_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_POTENTIAL_FACTORY_H_
