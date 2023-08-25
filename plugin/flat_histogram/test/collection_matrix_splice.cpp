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
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/energy.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wltm.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/collection_matrix_splice.h"

namespace feasst {

MonteCarlo monte_carlo2(const int thread, const int min, const int max,
    const int soft_min, const int soft_max) {
  DEBUG("min " << min);
  DEBUG("max " << max);
  const int trials_per = 1e2;
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  //mc.set(MakeRandomMT19937({{"seed", "1646430438"}}));
  mc.add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.fstprt"},
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
    //MakeWLTM({{"min_sweeps", "100000"}, {"new_sweep", "1"}, {"min_flatness", "25"}, {"collect_flatness", "20"}})));//, {"max_block_operations", "6"}})));
    MakeTransitionMatrix({{"min_sweeps", "100000"}, {"new_sweep", "1"}})));//, {"max_block_operations", "6"}})));
  mc.add(MakeCheckEnergy({{"trials_per_write", str(trials_per)}}));
  mc.add(MakeTune({{"trials_per_write", str(trials_per)}, {"multistate", "true"}, {"file_name", "tune" + str(thread)}}));
  mc.add(MakeLogAndMovie({{"trials_per_write", str(trials_per)},
    {"file_name", "tmp/clones" + str(thread)}}));
  mc.add(MakeCriteriaUpdater({{"trials_per_update", str(trials_per)}}));
  mc.add(MakeCriteriaWriter({
    {"trials_per_write", str(trials_per)},
    {"file_name", "tmp/clones" + str(thread) + "_crit.txt"}}));
  mc.set(MakeCheckpoint({{"num_hours", "0.0001"},
    {"file_name", "tmp/clone" + str(thread) + ".fst"}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/clone_energy" + str(thread)},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}}));
  return mc;
}

CollectionMatrixSplice make_splice(const int max, const int min = 0) {
  auto cm = MakeCollectionMatrixSplice({{"min_window_size", "2"},
    {"ln_prob_file", "tmp/lnpi.txt"},
    {"ln_prob_file_append", "true"},
    {"hours_per", "0.00001"}});
  std::vector<std::vector<int> > bounds = WindowExponential({
    {"maximum", str(max)},
    {"minimum", str(min)},
    {"num", "2"},
    {"overlap", "0"},
    {"alpha", "0.7"}}).boundaries();
  DEBUG(feasst_str(bounds));
  for (int index = 0; index < static_cast<int>(bounds.size()); ++index) {
    const std::vector<int> bound = bounds[index];
    DEBUG(bound[0] << " " << bound[1]);
    auto clone = std::make_shared<MonteCarlo>(monte_carlo2(index, min, max, bound[0], bound[1]));
    cm->add(clone);
  }
  return test_serialize(*cm);
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
  CollectionMatrix cm = clones2.collection_matrix(0);
  LnProbability ln_prob = clones2.ln_prob();
  //DEBUG(feasst_str(clones2.collection_matrix(0).matrix()));
  //DEBUG(feasst_str(clones2.collection_matrix(1).matrix()));
  //DEBUG(feasst_str(clones2.collection_matrix().matrix()));
  DEBUG("macro " << clones2.flat_histogram(0).macrostate().soft_max());
  DEBUG("macro " << clones2.flat_histogram(1).macrostate().soft_min());
  //DEBUG("visits " << clones2.flat_histogram(0).bias().visits(5));
  //DEBUG("visits " << clones2.flat_histogram(1).bias().visits(5));
  clones2.adjust_bounds();
  DEBUG("macro " << clones2.flat_histogram(0).macrostate().soft_max());
  DEBUG("macro " << clones2.flat_histogram(1).macrostate().soft_min());
  //DEBUG("visits " << clones2.flat_histogram(0).bias().visits(5));
  //DEBUG("visits " << clones2.flat_histogram(1).bias().visits(5));
  //DEBUG(feasst_str(clones2.collection_matrix(0).matrix()));
  //DEBUG(feasst_str(clones2.collection_matrix(1).matrix()));
  //DEBUG(feasst_str(clones2.collection_matrix().matrix()));
  clones2.write("tmp/ln_prob.txt");
}

TEST(CollectionMatrixSplice, lj_fh_LONG) {
  CollectionMatrixSplice clones2 = make_splice(5, 1);
  clones2.get_clone(0)->write_to_file();
  while (!clones2.are_all_complete()) {
    clones2.run(0.0001);
    DEBUG("swap");
    clones2.adjust_bounds();
  }
  clones2.write("tmp/ln_prob.txt");
  clones2.get_clone(0)->write_to_file();
//  clones2.run_until_all_are_complete();
  LnProbability lnpi = clones2.ln_prob();
  EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.04);
  EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.04);
  EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.04);
  EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.025);
  EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.01);
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
  EXPECT_NEAR(clones2.clone(0).analyze(en_index[0]).analyze(1).accumulator().average(), -0.030574223333333334, 0.001*3);
  EXPECT_NEAR(clones2.clone(0).analyze(en_index[0]).analyze(2).accumulator().average(), -0.089928316, 0.002*3);
  EXPECT_NEAR(clones2.clone(1).analyze(en_index[0]).analyze(3).accumulator().average(), -0.1784570533333333, 0.004*3);
  EXPECT_NEAR(clones2.clone(1).analyze(en_index[0]).analyze(4).accumulator().average(), -0.29619201333333334, 0.006*3);
//  INFO(SeekAnalyze().reference("Energy", clones2.clone(1)).accumulator().average());

  MakeCheckpoint({{"file_name", "tmp/clones.fst"}})->write(clones2);
  auto clones3 = MakeCollectionMatrixSplice("tmp/clones.fst");
}

}  // namespace feasst
