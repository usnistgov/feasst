#include <cmath>
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/lennard_jones.h"
#include "system/include/hard_sphere.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_bond.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/visit_model_cell.h"
#include "system/include/model_two_body_table.h"
#include "monte_carlo/include/metropolis.h"
#include "charge/include/utils.h"
#include "charge/include/ewald.h"
#include "charge/include/charge_screened.h"
#include "charge/include/charge_screened_intra.h"
#include "charge/include/charge_self.h"
//#include "charge/include/slab_correction.h"
#include "charge/include/trial_add_multiple.h"
#include "charge/include/trial_remove_multiple.h"

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
