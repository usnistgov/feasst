
#ifndef FEASST_CORE_PHYSICAL_CONSTANTS_H_
#define FEASST_CORE_PHYSICAL_CONSTANTS_H_

namespace feasst {

/**
 * constants from CODATA
 * http://dx.doi.org/10.1103/RevModPhys.84.1527
 * http://dx.doi.org/10.1063/1.4724320
 */
// HWH be more precise about the units that are used
// and the ability to change the units later.
// HWH consider the new SI physical constants.

/// Boltzman constant in units of Joules per Kelvin
constexpr double boltzmann_constant = 1.3806488E-23;

/// Avogadro's constant units of number of particles per mol
constexpr double avogadro_constant = 6.02214129E+23;

/// Ideal gas constant in units of Joules per Kelvin per mol
constexpr double ideal_gas_constant = boltzmann_constant*avogadro_constant;

/// Permitivity of vacuum in units of C^2/J/m
constexpr double permitivity_vacuum = 8.854187817E-12;

/// Elementary charge in units of C
constexpr double elementary_charge = 1.602176565E-19;

}  // namespace feasst

#endif  // FEASST_CORE_PHYSICAL_CONSTANTS_H_
