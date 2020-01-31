#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "opt_lj/include/visit_model_opt_lj.h"
#include "system/include/energy_map_all.h"

namespace feasst {

TEST(MonteCarlo, NVT_opt_lj_benchmark) {
  MonteCarlo mc;
  mc_lj(&mc);
  // HWH Perhaps implement as optimized potential instead of optimized VisitModel
  mc.add_to_optimized(Potential(MakeLennardJones(), //HWH: prevents ModelEmpty... how to remove?
                                MakeVisitModelOptLJ()));
 //                               MakeVisitModelOptLJ(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.seek_num_particles(50);
  // mc.seek_num_particles(250);
  // mc.attempt(1e6);  // 4.1 sec with 50 [5 sec 11/21/19, 4.3 on 1/6/19, 4.5 on 1/30/19 (4.18 VERBOSE_LEVEL=0)]
  // mc.seek_num_particles(450);
  // mc.attempt(1e5);  // 15 sec with 450 on slow computer
  mc.attempt(1e3);
  // DEBUG("\n" << mc.timer_str());
}

}
