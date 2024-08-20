
#ifndef FEASST_CHARGE_UTILS_H_
#define FEASST_CHARGE_UTILS_H_

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

}  // namespace feasst

#endif  // FEASST_CHARGE_UTILS_H_
