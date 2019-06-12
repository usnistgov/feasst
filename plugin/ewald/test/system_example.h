
#ifndef FEASST_EWALD_SYSTEM_EXAMPLE_H_
#define FEASST_EWALD_SYSTEM_EXAMPLE_H_

#include "system/include/system.h"
#include "ewald/include/utils_ewald.h"
#include "system/include/model_lj.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

inline System spce() {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  System sys;
  add_ewald_with(MakeModelLJ(), &config, &sys);
  sys.add(config);
  { Potential potential;
    potential.set_visit_model(MakeLongRangeCorrections());
    sys.add(potential);
  }
  sys.precompute();
  return sys;
}

}  // namespace feasst

#endif  // FEASST_EWALD_SYSTEM_EXAMPLE_H_
