
#ifndef FEASST_EWALD_UTILS_H_
#define FEASST_EWALD_UTILS_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

/**
Convert a temperature in units of Kelvin to energy in units of kJ/mol
after multiplying by the ideal gas constant.
Use the default physical constants (as given in ModelParams).
*/
double kelvin2kJpermol(const double kelvin);

/// Same as above, but obtain the physical constants from a Configuration.
double kelvin2kJpermol(const double kelvin, const Configuration& config);

/**
Prepare a MonteCarlo for a typical SPC/E simulation.

args:
- cubic_box_length: default 20
- particle: default forcefield/data.spce
- alphaL: Ewald alpha parameter (default: 5.6)
- kmax_squared: maximum squared wave vector (default: 38)
- physical_constants: as described in Configuration (default: CODATA2018)
- xyz_file: optionally load xyz file.
- dual_cut: add cell list with this width and also create a short-range
  reference potential with this cutoff.
  Ignore if -1 (default: -1).
 */
System spce(const argtype& args = argtype());

/**
Prepare a MonteCarlo for a typical RPM simulation.

args:
- cubic_box_length: default 12
- particle0: default plugin/ewald/forcefield/data.rpm_plus
- particle1: default plugin/ewald/forcefield/data.rpm_minus
- cutoff: (default: None, use particle file)
- alphaL: Ewald alpha parameter (default: 5.6)
- kmax_squared: maximum squared wave vector (default: 38)
- beta_mu: dimensionless (default: -13.94)
- delta: optional size disparity, sigma_+/- = 1 +/- delta
- charge_ratio: ratio of positive to negative charge (default: 1)
- dual_cut: add cell list with this width and also create a short-range
  reference potential with this cutoff.
  Ignore if -1 (default: -1).
 */
System rpm(const argtype& args = argtype());

}  // namespace feasst

#endif  // FEASST_EWALD_UTILS_H_
