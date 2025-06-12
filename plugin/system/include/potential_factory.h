
#ifndef FEASST_SYSTEM_POTENTIAL_FACTORY_H_
#define FEASST_SYSTEM_POTENTIAL_FACTORY_H_

#include <memory>
#include <vector>
#include <string>
// #include "utils/include/timer.h"

namespace feasst {

class Configuration;
class Potential;
class Select;

/**
  A list of potentials.

  By default, if one of the potentials returns a large energy, then the rest
  of the potentials are not calculated.
  This is the default behavior, which can be disabled with remove_opt_overlap()
  as done in MayerSampling.
 */
class PotentialFactory {
 public:
  PotentialFactory() {}

  /// Return the user name.
  const std::string& user_name() const { return user_name_; }

  /// Set the user name.
  void set_user_name(const std::string& name) { user_name_ = name; }

  /// Add potential.
  void add(std::shared_ptr<Potential> potential);

  /// Set a potential.
  void set(const int index, std::shared_ptr<Potential> potential);

  /// Return the potentials.
  const std::vector<std::shared_ptr<Potential> >& potentials() const {
    return potentials_; }

  /// Return a potential by index of order added.
  const Potential& potential(const int index) const {
    return const_cast<Potential&>(*potentials_[index]); }

  // Return a potential by index of order added.
  Potential * get_potential(const int index) {
    return potentials_[index].get(); }

  /// Return the number of potentials.
  int num() const { return static_cast<int>(potentials_.size()); }

  /// Precompute quantities for optimizations before calculation of energies.
  void precompute(Configuration * config);

  /// Precompute a particular potential by index.
  void precompute(const int index, Configuration * config);

  /// Remove optimization when overlap is detected, which is default.
  void remove_opt_overlap();

  /// Compute the energy of the given configuration.
  double energy(Configuration * config);

  /// Compute the energy of the selection in the configuration.
  double select_energy(const Select& select, Configuration * config);

  /// Return the profile of energies that were last computed.
  std::vector<double> stored_energy_profile() const;

  /// Return the last computed value of the energy.
  double stored_energy() const;

  /// Return a human-readable status.
  std::string str() const;

  /// Change the volume.
  void change_volume(const double delta_volume, const int dimension,
    Configuration * config);

  /// Revert any changes to the configuration due to the last energy computation
  void revert(const Select& select);

  /// Finalize changes due to perturbations.
  void finalize(const Select& select, Configuration * config);

  /// Set cache to load energy calculations.
  void load_cache(const bool load);

  /// Set cache to unload energy calclatuions.
  void unload_cache(const PotentialFactory& factory);

  void synchronize_(const PotentialFactory& factory, const Select& perturbed);

  void check(const Configuration& config) const;

  /// Serialize.
  void serialize(std::ostream& sstr) const;

  /// Deserialize.
  explicit PotentialFactory(std::istream& sstr);
  virtual ~PotentialFactory() {}

 private:
  std::vector<std::shared_ptr<Potential> > potentials_;
  int opt_overlap_ = 1;
//  Timer timer_;
  std::string user_name_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_POTENTIAL_FACTORY_H_
