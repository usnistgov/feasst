#include "configuration/include/configuration.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/model_params.h"
#include "charge/include/utils.h"

namespace feasst {

double kelvin2kJpermol(const double kelvin) {
  ModelParams model_params;
  const double R = model_params.physical_constants().ideal_gas_constant();
  return kelvin*R/1000.;
}

double kelvin2kJpermol(const double kelvin, const Configuration& config) {
  const double R = config.physical_constants().ideal_gas_constant();
  return kelvin*R/1000.;
}

}  // namespace feasst
