#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "prefetch/include/prefetch.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/utils.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/flat_histogram.h"

namespace feasst {

TEST(Prefetch, NVT_benchmark) {
  Prefetch mc;
  lennard_jones(&mc, {{"steps_per", str(1e1)}});
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  // mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.activate_prefetch(false);
  SeekNumParticles(50).with_trial_add().run(&mc);
  // activate prefetch after initial configuration
  mc.activate_prefetch(true);
  // mc.attempt(1e6);  // ~3.5 seconds (now 4.1)
  mc.attempt(1e2);
  EXPECT_EQ(mc.analyze(0).steps_since_write(),
            mc.modify(0).steps_since_update());
}

TEST(Prefetch, MUVT) {
  auto mc = MakePrefetch({{"steps_per_check", "1"}});
  lennard_jones(mc.get(), {{"steps_per", str(1e1)}});
  add_common_steppers(mc.get(), {{"steps_per", str(1e4)},
                                 {"file_append", "tmp/lj"}});
  //mc_lj(mc.get(), 8, "../forcefield/data.lj", 1e1, true, false);
  // mc->set(MakeRandomMT19937({{"seed", "default"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578665877"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578667496"}}));
  mc->set(MakeRandomMT19937({{"seed", "1804289383"}}));
  add_trial_transfer(mc.get(), {{"particle_type", "0"}});
  // mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-2"}}));
  mc->set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}})),
    MakeTransitionMatrix({{"min_sweeps", "10"}}),
    {{"beta", str(1./1.5)},
     {"chemical_potential", "-2.352321"}}));
  mc->add(MakeCriteriaUpdater({{"steps_per", str(1e1)}}));

//  // initialize ghosts the same
//  mc->seek_num_particles(100);
//  mc->get_system()->get_configuration()->remove_particles(
//    mc->configuration().selection_of_all());

  mc->activate_prefetch(true);
  mc->attempt(1);
  std::vector<std::shared_ptr<FlatHistogram> > fhs(mc->pool().size());
  std::vector<std::shared_ptr<TransitionMatrix> > tms(mc->pool().size());
  for (int trial = 0; trial < 1e2; ++trial) {
    mc->attempt(1);
    // INFO(mc->pool().size());
    for (int ipool = 0; ipool < static_cast<int>(mc->pool().size()); ++ipool) {
      std::stringstream ss;
      mc->clone_(ipool)->criteria().serialize(ss);
      fhs[ipool] = std::make_shared<FlatHistogram>(ss);
      std::stringstream ss2;
      fhs[ipool]->bias().serialize(ss2);
      tms[ipool] = std::make_shared<TransitionMatrix>(ss2);
      if (ipool != 0) {
        ASSERT(fhs[0]->is_equal(*fhs[ipool], NEAR_ZERO), "hi");
        ASSERT(tms[0]->is_equal(*tms[ipool], 1e-8), "hi");
      }
    }
  }

  Prefetch mc2 = test_serialize(*mc);
  EXPECT_EQ(1, mc2.steps_per_check());
}

}  // namespace feasst
