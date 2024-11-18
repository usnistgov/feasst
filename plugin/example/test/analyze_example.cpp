#include <cmath>
#include "utils/test/utils.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/monte_carlo.h"
#include "example/include/analyze_example.h"

namespace feasst {

TEST(AnalyzeExample, serialize) {
  auto obj = MakeAnalyzeExample({{"output_file", "tmp"}});
  auto obj2 = test_serialize_unique(*obj);
}

TEST(AnalyzeExample, ideal_gas_fluid_geometric_center_LONG) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../particle/atom.fstprt"},
                       {"cubic_side_length", "8"},
                       {"add_particles_of_type0", "100"}}},
    {"Potential", {{"Model", "IdealGas"}}},
    {"ThermoParams", {{"beta", "1"}}},
    {"Metropolis", {{"num_trials_per_iteration", "1e2"}}},
    {"TrialTranslate", {{"tunable_param", "4."}}},
    {"AnalyzeExample", {{"trials_per_write", str(1e3)},
                        {"output_file", "tmp/ig_center.csv"},
                        {"start_after_iteration", "1"}}},
    {"Run", {{"num_trials", "1e5"}}},
  }}, true);
  std::stringstream ss;
  mc->analyze(0).serialize(ss);
  AnalyzeExample analyze_example(ss);
  EXPECT_EQ(analyze_example.class_name(), "AnalyzeExample");
  for (int dim = 0; dim < mc->configuration().dimension(); ++dim) {
    const Accumulator& center = analyze_example.geometric_center(dim);
    EXPECT_TRUE(std::abs(center.average()) < 3*center.block_stdev());
  }
}

}  // namespace feasst
