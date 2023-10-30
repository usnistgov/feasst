
#ifndef FEASST_SYSTEM_THERM_PARAMS_H_
#define FEASST_SYSTEM_THERM_PARAMS_H_

#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include "utils/include/arguments.h"

namespace feasst {

/**
  Intensive thermodynamic state variables.
  Units should be consistent with the
  <a href="../../../particle/README.html">particle</a>.
 */
class ThermoParams {
 public:
  //@{
  /** @name Arguments
    - beta: inverse temperature, \f$ \beta = \frac{1}{k_B T} \f$.
    - pH: negative of the log-base-10 of the proton concentration.
    - chemical_potential[i]: chemical potential of the i-th particle type.
      The [i] is to be substituted for an integer 0, 1, 2, ...
      If only one particle type, you can drop the [i].
      The chemical potential must have the inverse units of \f$\beta\f$.
    - pressure: imposed isotropic system pressure.
   */
  explicit ThermoParams(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  /// Set beta.
  void set_beta(const double beta);

  /// Return beta.
  double beta() const;

  /// Set the pH.
  void set_pH(const double pH);

  /// Return the pH.
  double pH() const;

  /// Add a chemical potential for a given type of particle.
  /// Note that z has units length^{-dimension} such that Vz/N is unitless.
  // HWH Note: consider interfacing somehow with system/particle, etc
  // Perhaps System should contain criteria or MC kernel new object
  void add_chemical_potential(const double chemical_potential) {
    chemical_potentials_.push_back(chemical_potential); }

  /// Set the chemical potential of a given type.
  void set_chemical_potential(const double mu, const int particle_type = 0) {
    chemical_potentials_[particle_type] = mu; }

  /// Return the chemical potential of the particle type.
  double chemical_potential(const int particle_type = 0) const;

  /// Return the dimensionless product of beta and the chemical potential.
  double beta_mu(const int particle_type = 0) const;

  /// Return the pressure
  double pressure() const;

  /// Set the pressure.
  void set_pressure(const double pressure);

  /// Return a human readable string.
  std::string str() const;

  /// Return true if equivalent.
  bool is_equal(const ThermoParams& thermo_params) const;

  void serialize(std::ostream& ostr) const;
  explicit ThermoParams(std::istream& istr);

  //@}
 private:
  double beta_ = 0.;
  bool beta_initialized_ = false;
  std::vector<double> chemical_potentials_;
  double pH_ = 0.;
  bool pH_initialized_ = false;
  bool pressure_initialized_ = false;
  double pressure_ = 0.;
};

inline std::shared_ptr<ThermoParams> MakeThermoParams(
    argtype args = argtype()) {
  return std::make_shared<ThermoParams>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_THERM_PARAMS_H_
