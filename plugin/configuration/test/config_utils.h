
#ifndef FEASST_TEST_CONFIGURATION_UTILS_H_
#define FEASST_TEST_CONFIGURATION_UTILS_H_

#include "utils/include/arguments.h"
#include "configuration/include/configuration.h"

namespace feasst {

/// Return a LJ configuration from the SRSW sample 4.
Configuration lj_sample4();

/// Return an SPCE configuration from SRSW sample 1.
Configuration spce_sample1();

/**
  Return a configuration with two single-site particles.

  args:
  - cubic_box_length: default 6
 */
Configuration two_particle_configuration(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_TEST_CONFIGURATION_UTILS_H_
