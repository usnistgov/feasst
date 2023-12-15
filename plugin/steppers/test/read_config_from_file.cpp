#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/physical_constants.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/read_config_from_file.h"

namespace feasst {

TEST(ReadConfigFromFile, test) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {
      {"particle_type0", "../plugin/steppers/test/data/mab.fstprt"},
      {"xyz_euler_file", "../plugin/steppers/test/data/nvt0.xyze"},
    }},
    {"Potential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "1"}}},
    {"Metropolis", {{}}},
    {"ReadConfigFromFile", {{"input_file", "../plugin/steppers/test/data/nvt0.xyze"}, {"euler", "true"}}},
  }});
  const Configuration& config = mc->configuration();
  EXPECT_NEAR(-101.14905, config.particle(0).site(0).position().coord(0), NEAR_ZERO);
  EXPECT_NEAR(0.96721814, config.particle(0).site(0).euler().phi(), NEAR_ZERO);
  mc->attempt(1);
  EXPECT_FALSE(mc->criteria().is_complete());
  EXPECT_NEAR(-9.6221562, config.particle(0).site(0).position().coord(0), NEAR_ZERO);
  EXPECT_NEAR(1.2059475, config.particle(0).site(0).euler().theta(), NEAR_ZERO);
  mc->attempt(1);
  EXPECT_FALSE(mc->criteria().is_complete());
  EXPECT_NEAR(215.82084, config.particle(0).site(0).position().coord(0), NEAR_ZERO);
  EXPECT_NEAR(2.6758224, config.particle(0).site(4).euler().psi(), NEAR_ZERO);
  mc->attempt(1);
  EXPECT_TRUE(mc->criteria().is_complete());
}

}  // namespace feasst
