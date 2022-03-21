#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/energy.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/collection_matrix_splice.h"

namespace feasst {

MonteCarlo monte_carlo2(const int thread, const int min, const int max,
    const int soft_min, const int soft_max) {
  DEBUG("min " << min);
  DEBUG("max " << max);
  const int steps_per = 1e2;
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1646430438"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "8"},
                            {"particle_type0", "../forcefield/lj.fstprt"},
                            {"add_particles_of_type0", "1"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
  mc.run(MakeRun({{"until_num_particles", str(soft_min)}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}}),
      {{"soft_macro_min", str(soft_min)}, {"soft_macro_max", str(soft_max)}}),
    MakeTransitionMatrix({{"min_sweeps", "100000"}, {"new_sweep", "1"}})));//, {"max_block_operations", "6"}})));
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

CollectionMatrixSplice make_splice(const int max, const int min = 0) {
  CollectionMatrixSplice cmsplice;
  std::vector<std::vector<int> > bounds = WindowExponential({
    {"maximum", str(max)},
    {"minimum", str(min)},
    {"num", "2"},
    {"overlap", "0"},
    {"alpha", "0.7"}}).boundaries();
  for (int index = 0; index < static_cast<int>(bounds.size()); ++index) {
    const std::vector<int> bound = bounds[index];
    DEBUG(bound[0] << " " << bound[1]);
    auto clone = std::make_shared<MonteCarlo>(monte_carlo2(index, min, max, bound[0], bound[1]));
    cmsplice.add(clone);
  }
  return test_serialize(cmsplice);
}

TEST(CollectionMatrixSplice, lj_fh) {
  CollectionMatrixSplice clones2 = make_splice(12);
  EXPECT_EQ(clones2.num(), 2);
  DEBUG("num " << clones2.clone(0).configuration().num_particles());
  DEBUG("num " << clones2.clone(1).configuration().num_particles());
  DEBUG("macro " << clones2.flat_histogram(0).macrostate().histogram().str());
  DEBUG("macro " << clones2.flat_histogram(1).macrostate().histogram().str());
  EXPECT_FALSE(clones2.are_all_complete());
  clones2.run(0.0001);
  EXPECT_FALSE(clones2.are_all_complete());
  TripleBandedCollectionMatrix cm = clones2.collection_matrix(0);
  LnProbability ln_prob = clones2.ln_prob();
  DEBUG(feasst_str(clones2.collection_matrix(0).matrix()));
  DEBUG(feasst_str(clones2.collection_matrix(1).matrix()));
  DEBUG(feasst_str(clones2.collection_matrix().matrix()));
  DEBUG("macro " << clones2.flat_histogram(0).macrostate().soft_max());
  DEBUG("macro " << clones2.flat_histogram(1).macrostate().soft_min());
  DEBUG("visits " << clones2.flat_histogram(0).bias().visits(5));
  DEBUG("visits " << clones2.flat_histogram(1).bias().visits(5));
  clones2.adjust_bounds(1);
  DEBUG("macro " << clones2.flat_histogram(0).macrostate().soft_max());
  DEBUG("macro " << clones2.flat_histogram(1).macrostate().soft_min());
  DEBUG("visits " << clones2.flat_histogram(0).bias().visits(5));
  DEBUG("visits " << clones2.flat_histogram(1).bias().visits(5));
  DEBUG(feasst_str(clones2.collection_matrix(0).matrix()));
  DEBUG(feasst_str(clones2.collection_matrix(1).matrix()));
  DEBUG(feasst_str(clones2.collection_matrix().matrix()));
}

TEST(CollectionMatrixSplice, lj_fh_LONG) {
  CollectionMatrixSplice clones2 = make_splice(5, 1);
  while (!clones2.are_all_complete()) {
    clones2.run(0.0001);
    DEBUG("swap");
    clones2.adjust_bounds(1);
  }
  LnProbability lnpi = clones2.ln_prob();
  EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.02);
  EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.02);
  EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.02);
  EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.01);
  EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.005);
  std::vector<double> en0 = SeekAnalyze().multistate_data("Energy", clones2.clone(0));
  DEBUG(feasst_str(en0));
  std::vector<double> en1 = SeekAnalyze().multistate_data("Energy", clones2.clone(1));
  DEBUG(feasst_str(en1));
  const std::vector<int> en_index = SeekAnalyze().index("Energy", clones2.clone(0));
  DEBUG(feasst_str(en_index));
//  INFO(clones2.clone(0).analyze(en_index[0]).analyze(0).accumulator().average());
//  INFO(clones2.clone(0).analyze(en_index[0]).analyze(1).accumulator().average());
//  INFO(clones2.clone(0).analyze(en_index[0]).analyze(2).accumulator().average());
//  INFO(clones2.clone(1).analyze(en_index[0]).analyze(2).accumulator().average());
//  INFO(clones2.clone(1).analyze(en_index[0]).analyze(3).accumulator().average());
//  INFO(clones2.clone(1).analyze(en_index[0]).analyze(4).accumulator().average());
  EXPECT_NEAR(clones2.clone(0).analyze(en_index[0]).analyze(0).accumulator().average(), -0.000605740233333333, 1e-8);
  EXPECT_NEAR(clones2.clone(0).analyze(en_index[0]).analyze(1).accumulator().average(), -0.030574223333333334, 0.001*2);
  EXPECT_NEAR(clones2.clone(0).analyze(en_index[0]).analyze(2).accumulator().average(), -0.089928316, 0.002*2);
  EXPECT_NEAR(clones2.clone(1).analyze(en_index[0]).analyze(3).accumulator().average(), -0.1784570533333333, 0.004*2);
  EXPECT_NEAR(clones2.clone(1).analyze(en_index[0]).analyze(4).accumulator().average(), -0.29619201333333334, 0.006*2);
//  INFO(SeekAnalyze().reference("Energy", clones2.clone(1)).accumulator().average());
}

}  // namespace feasst
