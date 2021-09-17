
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
  Return a system for a typical SPC/E simulation.

  args:
  - Configuration args (defaults; cubic_box_length: 20,
    particle_type: /path/to/feasst/forcefield/spce.fstprt)
  - lrc: use long range corrections (default: true)
  - alphaL: Ewald alpha parameter (default: 5.6)
  - kmax_squared: maximum squared wave vector (default: 38)
  - physical_constants: as described in Configuration (default: CODATA2018)
  - xyz_file: optionally load xyz file.
  - dual_cut: add cell list with this width and also create a short-range
    reference potential with this cutoff.
    Ignore if -1 (default: -1).
  - table_size: size of tabular potential to speed up erfc (default: 1e6).
 */
System spce(argtype args = argtype());

/**
Prepare a MonteCarlo for a typical RPM simulation.

args:
- cubic_box_length: default 12
- particle0: default plugin/charge/forcefield/rpm_plus.fstprt
- particle1: default plugin/charge/forcefield/rpm_minus.fstprt
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
System rpm(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_EWALD_UTILS_H_
