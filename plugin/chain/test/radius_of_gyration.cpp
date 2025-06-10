#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/physical_constants.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/read_config_from_file.h"
#include "chain/include/radius_of_gyration.h"

namespace feasst {

TEST(RadiusOfGyration, test) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {
      {"particle_type0", "../plugin/steppers/test/data/mab.txt"},
      {"xyz_euler_file", "../plugin/steppers/test/data/nvt0.xyze"},
    }},
    {"Potential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1"}}},
    {"Metropolis", {{}}},
    {"RadiusOfGyration", {{"print_histogram", "true"}, {"min", "27"}, {"max", "27"}, {"width", "1"}, {"output_file", "tmp/rg.csv"}}},
    {"ReadConfigFromFile", {{"input_file", "../plugin/steppers/test/data/nvt0.xyze"}, {"euler", "true"}}},
    {"Run", {{"until", "complete"}}},
    {"WriteStepper", {{"analyze_name", "RadiusOfGyration"}}},
  }}, true);
  EXPECT_NEAR(27.336812419903538, mc->analyze(0).accumulator().average(), 1e-8);
}

}  // namespace feasst
