#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "math/include/histogram.h"
#include "math/include/random_mt19937.h"
#include "math/include/histogram.h"
#include "configuration/include/configuration.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/clones.h"

namespace feasst {

std::unique_ptr<MonteCarlo> monte_carlo(const int thread, const int min, const int max) {
  const int trials_per = 1e2;
  auto mc = std::make_unique<MonteCarlo>();
  //mc->set(MakeRandomMT19937({{"seed", "1635444301"}}));
  mc->add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.fstprt"},
                            {"add_particles_of_type0", "1"}}));
  mc->add(MakePotential(MakeLennardJones()));
  mc->add(MakePotential(MakeLongRangeCorrections()));
  mc->set(MakeThermoParams({{"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
  mc->run(MakeRun({{"until_num_particles", str(min)}}));
  mc->set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "10"}})));//, {"max_block_operations", "6"}})));
  mc->add(MakeCheckEnergy({{"trials_per_update", str(trials_per)}}));
  mc->add(MakeTune());
//  mc->add(MakeLogAndMovie({{"trials_per_write", str(trials_per)},
//    {"output_file", "tmp/clones" + str(thread)}}));
  mc->add(MakeCriteriaUpdater({{"trials_per_update", str(trials_per)}}));
  mc->add(MakeCriteriaWriter({
    {"trials_per_write", str(trials_per)},
    {"output_file", "tmp/clones" + str(thread) + "_crit.txt"}}));
  mc->set(MakeCheckpoint({{"num_hours", "0.0001"},
    {"checkpoint_file", "tmp/clone" + str(thread) + ".fst"}}));
  mc->add(MakeEnergy({
    {"output_file", "tmp/clone_energy" + str(thread)},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
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
Clones make_clones(const int max, const int min = 0, const int overlap = 4) {
  Clones clones;
  std::vector<std::vector<int> > bounds = WindowExponential({
    {"maximum", str(max)},
    {"minimum", str(min)},
    {"num", "2"},
    {"overlap", str(overlap)},
    {"alpha", "2"}}).boundaries();
  for (int index = 0; index < static_cast<int>(bounds.size()); ++index) {
    const std::vector<int> bound = bounds[index];
    DEBUG(bound[0] << " " << bound[1]);
    std::unique_ptr<MonteCarlo> mcu = monte_carlo(index, bound[0], bound[1]);
    std::shared_ptr<MonteCarlo> mcs = std::move(mcu);
//    clone->set(MakeRandomMT19937({{"seed", "123"}}));
    clones.add(mcs);
  }
  return test_serialize(clones);
}

TEST(Clones, lj_fh) {
  Clones clones2 = make_clones(12);
  EXPECT_EQ(clones2.num(), 2);
  DEBUG("num " << clones2.clone(0).configuration().num_particles());
  DEBUG("num " << clones2.clone(1).configuration().num_particles());
  clones2.initialize_and_run_until_complete({{"omp_batch", str(1e1)}});
  DEBUG("0: " << feasst_str(clones2.flat_histogram(0)->bias().ln_prob().values()));
  DEBUG("1: " << feasst_str(clones2.flat_histogram(1)->bias().ln_prob().values()));
  EXPECT_NEAR(clones2.ln_prob().value(0), -36.9, 0.7);
  MakeCheckpoint({{"checkpoint_file", "tmp/rstclone"}})->write(clones2);
  Clones clones3;
  MakeCheckpoint({{"checkpoint_file", "tmp/rstclone"}})->read(&clones3);
  EXPECT_TRUE(clones3.ln_prob().is_equal(clones2.ln_prob(), 1e-8));

  Histogram macrostates;
  std::vector<double> energy, energy0, energy1;
  clones3.stitch(&macrostates, &energy, "Energy");
  EXPECT_EQ(macrostates.center_of_bin(0), 0);
  EXPECT_EQ(macrostates.center_of_bin(10), 10);
  EXPECT_EQ(macrostates.center_of_bin(11), 11);
  EXPECT_EQ(macrostates.center_of_bin(12), 12);
  energy0 = SeekAnalyze().multistate_data("Energy", clones3.clone(0));
  energy1 = SeekAnalyze().multistate_data("Energy", clones3.clone(1));
  DEBUG(feasst_str(energy));
  DEBUG(feasst_str(energy0));
  DEBUG(feasst_str(energy1));
  for (int i = 0; i < 5; ++i) EXPECT_EQ(energy[i], energy0[i]);
  EXPECT_EQ(energy[5], 0.5*(energy0[5] + energy1[0]));
  EXPECT_EQ(energy[6], 0.5*(energy0[6] + energy1[1]));
  EXPECT_EQ(energy[7], 0.5*(energy0[7] + energy1[2]));
  EXPECT_EQ(energy[8], 0.5*(energy0[8] + energy1[3]));
  for (int i = 9; i < 13; ++i) EXPECT_EQ(energy[i], energy1[i - 5]);
}

double energy_av4(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

TEST(Clones, lj_fh_LONG) {
  Clones clones = make_clones(5, 1, 1);
  Clones clones2 = test_serialize(clones);
  clones2.initialize_and_run_until_complete(
    {{"omp_batch", str(1e5)}, {"ln_prob_file", "tmp/clones_fh.txt"}});
  for (int sweeps = 20; sweeps <= 1000; sweeps+=10) {
  //for (int sweeps = 20; sweeps <= 100; sweeps+=10) {
    DEBUG("sweeps: " << sweeps);
    auto clones3 = MakeClones("tmp/clone", 2, 0, ".fst");
    clones3->set_num_iterations_to_complete(sweeps);
    clones3->initialize_and_run_until_complete(
      {{"omp_batch", str(1e5)}, {"ln_prob_file", "tmp/clones_fh.txt"}});
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
  // removed extra overlap
  //EXPECT_NEAR(energy_av4(0, clones4->clone(1)), -0.000605740233333333, 1e-8);
  //EXPECT_NEAR(energy_av4(1, clones4->clone(1)), -0.030574223333333334, 0.001);
  //EXPECT_NEAR(energy_av4(2, clones4->clone(1)), -0.089928316, 0.002);
  EXPECT_NEAR(energy_av4(0, clones4->clone(1)), -0.1784570533333333, 0.0055);
  EXPECT_NEAR(energy_av4(1, clones4->clone(1)), -0.29619201333333334, 0.0075);
}

}  // namespace feasst
