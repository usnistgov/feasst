#include "utils/test/utils.h"
#include "monte_carlo/include/remove.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

TEST(Remove, serialize) {
  auto remove = std::make_shared<Remove>(argtype({{"name", "something"}}));
  Action remove2 = test_serialize(*remove);
}

TEST(Remove, comma_separated) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{}}},
    {"Log", {{}}},
    {"Movie", {{"output_file", "tmp/hi"}}},
  }}, true);
  EXPECT_EQ(2, mc->num_analyzers());
  MakeRemove({{"name", "Log,Movie"}})->run(mc.get());
  EXPECT_EQ(0, mc->num_analyzers());
}

}  // namespace feasst
