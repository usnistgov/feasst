
#ifndef FEASST_TEST_CHARGE_CHARGE_UTILS_H_
#define FEASST_TEST_CHARGE_CHARGE_UTILS_H_

#include "utils/include/arguments.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/model_params.h"
#include "system/include/system.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
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

// HWH deprecate
/*
  Return a system for a typical SPC/E simulation.

  args:
  - Configuration args with the following default values:
      - cubic_side_length: 20,
      - particle_type: /path/to/feasst/particle/spce.fstprt,
      - physical_constants: CODATA2018.
  - Ewald args.
  - lrc: use long range corrections (default: true)
  - dual_cut: add cell list with this width and also create a short-range
    reference potential with this cutoff.
    Ignore if -1 (default: -1).
  - erfc_table_size: size of tabular potential to speed up erfc (default: 2e4).
 */
System spce(argtype args = argtype());

// HWH deprecate
/*
  Return a system for a typical RPM simulation.

  args:
  - Configuration args with the following default values:
      - cubic_side_length: default 12
      - particle0: default plugin/charge/particle/rpm_plus.fstprt
      - particle1: default plugin/charge/particle/rpm_minus.fstprt
  - Ewald args.
  - cutoff: (default: None, use particle file)
  - delta: optional size disparity, sigma_+/- = 1 +/- delta
  - charge_ratio: ratio of positive to negative charge (default: 1)
  - dual_cut: add cell list with this width and also create a short-range
    reference potential with this cutoff.
    Ignore if -1 (default: -1).
 */
System rpm(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_TEST_CHARGE_CHARGE_UTILS_H_
