
#ifndef FEASST_EWALD_SYSTEM_EXAMPLE_H_
#define FEASST_EWALD_SYSTEM_EXAMPLE_H_

#include "system/include/system.h"
#include "ewald/include/utils_ewald.h"
#include "system/include/lennard_jones.h"
#include "system/include/potential.h"
#include "configuration/test/configuration_test.h"

namespace feasst {

inline System spce(const std::string physical_constants = "") {
  Configuration config = spce_sample();
  config.add_model_param("alpha", 5.6/config.domain().min_side_length());
  if (!physical_constants.empty()) {
    if (physical_constants == "CODATA2010") {
      config.set_physical_constants(MakeCODATA2010());
    } else if (physical_constants == "CODATA2018") {
      config.set_physical_constants(MakeCODATA2018());
    } else {
      ERROR("unrecognized");
    }
  }
  System sys;
  sys.add(config);
  add_ewald_with(MakeLennardJones(), &sys);
  sys.add(Potential(MakeLongRangeCorrections()));
  sys.precompute();
  return sys;
}

inline void test_cases(
    /// tuple of constants and expected energy values
    std::vector<std::tuple<std::shared_ptr<PhysicalConstants>, double> > cases,
    std::shared_ptr<Model> model,
    std::shared_ptr<VisitModel> visitor = MakeVisitModel()) {
  for (auto cse : cases) {
    INFO(std::get<0>(cse)->class_name());
    Configuration config = spce_sample();
    config.add_model_param("alpha", 5.6/config.domain().min_side_length());
    config.set_physical_constants(std::get<0>(cse));
    Potential potential(model, visitor);
    potential.precompute(&config);
    EXPECT_NEAR(std::get<1>(cse), potential.energy(&config), 1e-10);
    Potential potential2 = test_serialize(potential);
    EXPECT_NEAR(std::get<1>(cse), potential2.energy(&config), 1e-10);
  }
}

}  // namespace feasst

#endif  // FEASST_EWALD_SYSTEM_EXAMPLE_H_
