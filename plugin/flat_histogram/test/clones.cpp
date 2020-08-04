#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "math/include/random_mt19937.h"
#include "system/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_transfer.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/energy.h"
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
  mc.get_system()->get_configuration()->add_particle_of_type(0);
  mc.add(MakeFlatHistogram(
      MakeMacrostateNumParticles(
        Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
      MakeTransitionMatrix({{"min_sweeps", "10"}}),
      { {"beta", str(1./1.5)},
        {"chemical_potential", "-2.352321"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)},
      {"file_name", "tmp/clones" + str(thread)}}));
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/clones" + str(thread) + "_crit.txt"}}));
  mc.set(MakeCheckpoint({{"num_hours", "0.0001"},
    {"file_name", "tmp/clone" + str(thread) + ".fst"}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/clone_energy" + str(thread)},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}}));
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
Clones make_clones(const int max, const int min) {
  Clones clones;
  std::vector<std::vector<int> > bounds = WindowExponential({
    {"maximum", str(max)},
    {"minimum", str(min)},
    {"num", "2"},
    {"extra_overlap", "3"},
    {"alpha", "2"}}).boundaries();
  for (int index = 0; index < static_cast<int>(bounds.size()); ++index) {
    const std::vector<int> bound = bounds[index];
    INFO(bound[0] << " " << bound[1]);
    auto clone = std::make_shared<MonteCarlo>(monte_carlo(index, bound[0], bound[1]));
//    clone->set(MakeRandomMT19937({{"seed", "123"}}));
    if (index == 0 && bound[0] > 0) SeekNumParticles(bound[0]).run(clone.get());
    clones.add(clone);
  }
  return test_serialize(clones);
}

TEST(Clones, lj_fh) {
  Clones clones2 = make_clones(12, 0);
  EXPECT_EQ(clones2.num(), 2);
  DEBUG("num " << clones2.clone(0).configuration().num_particles());
  DEBUG("num " << clones2.clone(1).configuration().num_particles());
  clones2.initialize_and_run_until_complete({{"omp_batch", str(1e1)}});
  DEBUG("0: " << feasst_str(clones2.flat_histogram(0).bias().ln_prob().values()));
  DEBUG("1: " << feasst_str(clones2.flat_histogram(1).bias().ln_prob().values()));
  EXPECT_NEAR(clones2.ln_prob().value(0), -36.9, 0.7);
  MakeCheckpoint({{"file_name", "tmp/rstclone"}})->write(clones2);
  Clones clones3;
  MakeCheckpoint({{"file_name", "tmp/rstclone"}})->read(&clones3);
  EXPECT_TRUE(clones3.ln_prob().is_equal(clones2.ln_prob(), 1e-8));
}

double energy_av4(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

TEST(Clones, lj_fh_LONG) {
  Clones clones = make_clones(5, 1);
  Clones clones2 = test_serialize(clones);
  clones2.initialize_and_run_until_complete(
    {{"omp_batch", str(1e3)}, {"ln_prob_file", "tmp/clones_fh.txt"}});
  for (int sweeps = 20; sweeps <= 1000; sweeps+=10) {
    DEBUG("sweeps: " << sweeps);
    auto clones3 = MakeClones("tmp/clone", 2, 0, ".fst");
    clones3->set_num_iterations(sweeps);
    clones3->initialize_and_run_until_complete(
      {{"omp_batch", str(1e3)}, {"ln_prob_file", "tmp/clones_fh.txt"}});
  }

  auto clones4 = MakeClones("tmp/clone", 2, 0, ".fst");
  const LnProbability lnpi = clones4->ln_prob();
  EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.04);
  EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.04);
  EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.04);
  EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.04);
  EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.01);
  EXPECT_NEAR(energy_av4(0, clones4->clone(0)), -0.000605740233333333, 1e-8);
  EXPECT_NEAR(energy_av4(1, clones4->clone(0)), -0.030574223333333334, 0.001);
  EXPECT_NEAR(energy_av4(2, clones4->clone(0)), -0.089928316, 0.002);
  EXPECT_NEAR(energy_av4(3, clones4->clone(0)), -0.1784570533333333, 0.004);
  EXPECT_NEAR(energy_av4(0, clones4->clone(1)), -0.000605740233333333, 1e-8);
  EXPECT_NEAR(energy_av4(1, clones4->clone(1)), -0.030574223333333334, 0.001);
  EXPECT_NEAR(energy_av4(2, clones4->clone(1)), -0.089928316, 0.002);
  EXPECT_NEAR(energy_av4(3, clones4->clone(1)), -0.1784570533333333, 0.004);
  EXPECT_NEAR(energy_av4(4, clones4->clone(1)), -0.29619201333333334, 0.006);
}

}  // namespace feasst
