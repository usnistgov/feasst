#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "system/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_transfer.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/clones.h"

namespace feasst {

MonteCarlo monte_carlo(const int thread, const int min, const int max) {
  const int steps_per = 1e2;
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
  mc.add(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "10"}}),
    { {"beta", str(1./1.5)},
      {"chemical_potential", "-2.352321"}}));
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/clones" + str(thread) + "_crit.txt"},
  }));
  return mc;
}

// macrostates
// 0 1 2 3 4 5 6 7 8            : 9 total
//           5 6 7 8 9 10 11 12 : 6 total
// if TM (yes), ignore upper min and lower max
// 0 1 2 3 4 5 6 7           : 9 total
//             6 7 8 9 10 11 : 6 total
//
// for 3 windows
// 0 1 2 3 4 5 6 7              : 8 total
//         4 5 6 7 8 9 10       : 7 total
//               7 8 9 10 11 12 : 6 total
//
// now consider ignoring ends
// 0 1 2 3 4 5 6                : 8 total
//           5 6 7 8 9          : 7 total
//                 8 9 10 11 12 : 6 total
TEST(Clones, lj_fh) {
  Clones clones;
  std::vector<std::vector<int> > bounds = WindowExponential({
    {"maximum", "12"},
    {"num", "2"},
    {"extra_overlap", "3"},
    {"alpha", "2"}}).boundaries();
  for (int index = 0; index < static_cast<int>(bounds.size()); ++index) {
    const std::vector<int> bound = bounds[index];
    DEBUG(bound[0] << " " << bound[1]);
    auto clone = std::make_shared<MonteCarlo>(monte_carlo(index, bound[0], bound[1]));
    if (index == 0 && bound[0] > 0) SeekNumParticles(bound[0]).run(clone.get());
    clone->add(MakeLogAndMovie({{"steps_per", str(1e5)},
      {"file_name", "tmp/clones" + str(clones.num())}}));
    clone->add(MakeCheckEnergyAndTune({{"steps_per", str(1e5)}}));
    clones.add(clone);
  }
  Clones clones2 = test_serialize(clones);
  EXPECT_EQ(clones.num(), 2);
  //clones2.initialize();
  DEBUG("num " << clones2.clone(0).configuration().num_particles());
  DEBUG("num " << clones2.clone(1).configuration().num_particles());
  //clones2.run_until_complete();
  clones2.initialize_and_run_until_complete();
  DEBUG("0: " << feasst_str(clones2.flat_histogram(0).bias().ln_prob().values()));
  DEBUG("1: " << feasst_str(clones2.flat_histogram(1).bias().ln_prob().values()));
//  INFO(feasst_str(clones2.ln_prob().values()));
  EXPECT_NEAR(clones2.ln_prob().value(0), -36.9, 0.7);

  MakeCheckpoint({{"file_name", "tmp/rstclone"}})->write(clones2);
  Clones clones3;
  MakeCheckpoint({{"file_name", "tmp/rstclone"}})->read(&clones3);
//  INFO(feasst_str(clones3.ln_prob().values()));
  EXPECT_TRUE(clones3.ln_prob().is_equal(clones2.ln_prob(), 1e-8));
}

}  // namespace feasst
