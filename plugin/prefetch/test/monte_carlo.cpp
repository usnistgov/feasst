#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "prefetch/include/prefetch.h"
#include "flat_histogram/test/flat_histogram_test.h"
#include "steppers/include/criteria_updater.h"

namespace feasst {

TEST(Prefetch, NVT_benchmark) {
  Prefetch mc;
  mc_lj(&mc, 8, "../forcefield/data.lj", 1e1);
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.activate_prefetch(false);
  mc.seek_num_particles(50);
  // activate prefetch after initial configuration
  mc.activate_prefetch(true);
  // mc.attempt(1e6);  // ~3.5 seconds (now 4.1)
  mc.attempt(1e2);
  //INFO("num " << mc.trials().num_attempts());
}

TEST(Prefetch, MUVT) {
  auto mc = MakePrefetch({{"steps_per_check", "1"}});
  mc_lj(mc.get(), 8, "../forcefield/data.lj", 1e1);
  // mc->set(MakeRandomMT19937({{"seed", "default"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578665877"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578667496"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1579010457"}}));
  add_trial_transfer(mc.get(), {{"particle_type", "0"}});
  // mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-2"}}));
  mc->set(crit_fh(0));
  mc->add(MakeCriteriaUpdater({{"steps_per", str(1e1)}}));

//  // initialize ghosts the same
//  mc->seek_num_particles(100);
//  mc->get_system()->get_configuration()->remove_particles(
//    mc->configuration().selection_of_all());

  mc->activate_prefetch(true);
  mc->attempt(1e2);

  Prefetch mc2 = test_serialize(*mc);
  EXPECT_EQ(1, mc2.steps_per_check());
}

}  // namespace feasst
