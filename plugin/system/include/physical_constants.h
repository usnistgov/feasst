
#ifndef FEASST_SYSTEM_PHYSICAL_CONSTANTS_H_
#define FEASST_SYSTEM_PHYSICAL_CONSTANTS_H_

#include <math.h>
#include "math/include/constants.h"

namespace feasst {

/**
 * constants from CODATA
 *
 * 2014:
 * http://dx.doi.org/10.1103/RevModPhys.88.035009
 * http://dx.doi.org/10.1063/1.4954402
 *
 * 2010 (previous, but no longer in use in feasst):
 * http://dx.doi.org/10.1103/RevModPhys.84.1527
 * http://dx.doi.org/10.1063/1.4724320
 */

/// Boltzman constant in units of Joules per Kelvin
constexpr double boltzmann_constant = 1.3806491E-23; // CODATA 2014
//nstexpr double boltzmann_constant = 1.3806488E-23; // CODATA 2010

/// Avogadro's constant units of number of particles per mol
constexpr double avogadro_constant = 6.02214076E+23; // CODATA 2014
//nstexpr double avogadro_constant = 6.02214129E+23; // CODATA 2010

/// Ideal gas constant in units of Joules per Kelvin per mol
constexpr double ideal_gas_constant = boltzmann_constant*avogadro_constant;

/// Permitivity of vacuum in units of C^2/J/m
constexpr double permitivity_vacuum = 8.8541878128E-12;// CODATA 2014
//nstexpr double permitivity_vacuum = 8.854187817E-12; // CODATA 2010

/// Elementary charge in units of C
constexpr double elementary_charge = 1.602176634E-19; // CODATA 2014
//nstexpr double elementary_charge = 1.602176565E-19; // CODATA 2010

/// Convert e^2/Angstrom to kJ/mol by a factor of units (kJ*A/e^2/mol)
const double charge_conversion = pow(elementary_charge, 2)/
              (4*PI*permitivity_vacuum*1e3/1e10/avogadro_constant);

}  // namespace feasst

#endif  // FEASST_SYSTEM_PHYSICAL_CONSTANTS_H_
