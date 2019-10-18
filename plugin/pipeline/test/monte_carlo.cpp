#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "pipeline/include/pipeline.h"

namespace feasst {

TEST(Pipeline, NVT_benchmark) {
  Pipeline mc;
  mc_lj(&mc);
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  //mc_lj(&mc, 8, "../forcefield/data.lj", 2e0);
  mc.activate_pipeline(false);
  mc.seek_num_particles(50);
  //mc.seek_num_particles(250);
  //INFO("done");
  // mc.activate_pipeline(true);
  //mc.set(mc2.system());
  // mc.attempt(1e6);  // ~3.5 seconds (now 4.1)
  mc.attempt(1e2);
  //INFO("num " << mc.trials().num_attempts());
}

}  // namespace feasst
