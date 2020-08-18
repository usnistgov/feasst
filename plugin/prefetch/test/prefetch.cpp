#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "prefetch/include/prefetch.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/energy.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/check_properties.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/flat_histogram.h"
#include "ewald/include/utils.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/trial_avb2_half.h"
#include "cluster/include/trial_avb4.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/trial_transfer_avb.h"

namespace feasst {

void run_prefetch(const int trials, const int steps_per) {
  auto mc = MakePrefetch();
//  mc->set(MakeRandomMT19937({{"seed", "1592943710"}}));
  mc->set(MakeRandomMT19937({{"seed", "1596650884"}}));
  mc->set(lennard_jones());
  mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}, {"num_steps", "1"}}));
  mc->add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/lj"}}));
  mc->add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
  mc->activate_prefetch(false);
  SeekNumParticles(50).with_trial_add().run(mc.get());
  // activate prefetch after initial configuration
  mc->activate_prefetch(true);
  mc->attempt(trials);
  EXPECT_EQ(mc->analyze(0).steps_since_write(),
            mc->modify(0).steps_since_update());
}

TEST(Prefetch, NVT_benchmark) {
  run_prefetch(1e3, 1e1);
}

TEST(Prefetch, NVT_benchmark_LONG) {
  run_prefetch(1e6, 1e3); // 5.4s on 4 cores of i7-4770K @ 3.5GHz
}

TEST(Prefetch, MUVT) {
  auto mc = MakePrefetch({{"steps_per_check", "1"}});
  mc->set(lennard_jones());
  mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeLogAndMovie({{"steps_per", str(1e1)}, {"file_name", "tmp/lj"}}));
  mc->add(MakeCheckEnergyAndTune({{"steps_per", str(1e1)}}));
  //mc_lj(mc.get(), 8, "../forcefield/data.lj", 1e1, true, false);
  // mc->set(MakeRandomMT19937({{"seed", "default"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578665877"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578667496"}}));
  mc->set(MakeRandomMT19937({{"seed", "1804289383"}}));
  mc->add(MakeTrialAdd({{"particle_type", "0"}}));
  mc->add(MakeTrialRemove({{"particle_type", "0"}}));
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
        ASSERT(fhs[0]->is_equal(*fhs[ipool], 1e-8), "hi");
        ASSERT(tms[0]->is_equal(*tms[ipool], 1e-8), "hi");
      }
    }
  }

  Prefetch mc2 = test_serialize(*mc);
  EXPECT_EQ(1, mc2.steps_per_check());
}

TEST(Prefetch, NVT_spce) {
  auto mc = MakePrefetch({{"synchronize", "true"}});
  //auto mc = MakePrefetch({{"synchronize", "false"}});
  // mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->set(spce());
  const int steps_per = 1e2;
  mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/lj"}}));
  //mc->add(MakeCheckProperties({{"steps_per", "1"}}));
  mc->add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc->add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
  // mc->set(MakeRandomMT19937({{"seed", "default"}}));
  mc->activate_prefetch(false);
  SeekNumParticles(50).with_trial_add().run(mc.get());
  // activate prefetch after initial configuration
  mc->activate_prefetch(true);
  // mc->attempt(1e6);  // ~3.5 seconds (now 4.1)
  mc->attempt(1e2);
  EXPECT_EQ(mc->analyze(0).steps_since_write(),
            mc->modify(1).steps_since_update());
}

TEST(Prefetch, AVB) {
  const int steps_per = 1e2;
  auto monte_carlo = MakePrefetch({{"synchronize", "true"},
                                   {"steps_per_check", str(steps_per)}});
  monte_carlo->add(Configuration(MakeDomain({{"cubic_box_length", "6"}}),
                                {{"particle_type", "../forcefield/data.lj"}}));
  monte_carlo->add(Potential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  //monte_carlo->add(Potential(MakeLennardJones()));
  monte_carlo->add(MakeMetropolis({{"beta", "0.00001"}, {"chemical_potential", "50."}}));
  SeekNumParticles(50).with_trial_add().run(monte_carlo.get());
  monte_carlo->add(MakeMetropolis({{"beta", "0.2"}, {"chemical_potential", "-20."}}));
  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "3"},
                                                 {"minimum_distance", "1"}});
  // Something wrong with adding TrialFactories when using prefetch...
  // I think they don't end up selecting the right trial type
//  monte_carlo->add(MakeTrialAVB2(neighbor_criteria));
  monte_carlo->add(MakeTrialAVB2Half(neighbor_criteria, {{"out_to_in", "true"}}));
  monte_carlo->add(MakeTrialAVB2Half(neighbor_criteria, {{"out_to_in", "false"}}));
  monte_carlo->add(MakeTrialAVB4(neighbor_criteria));
//  monte_carlo->add(MakeTrialTransferAVB(neighbor_criteria,
//    {{"particle_type", "0"}}));
  monte_carlo->add(MakeTrialTranslate({{"tunable_param", "4"}}));
//  monte_carlo->add(MakeTrialTransfer({{"particle_type", "0"}}));
  monte_carlo->add(MakeLogAndMovie({{"steps_per", str(steps_per)},
                             {"file_name", "tmp/pljavb.xyz"},
                             {"clear_file", "true"}}));
  monte_carlo->add(MakeCheckEnergy({{"steps_per", str(steps_per)},
                                   {"tolerance", str(1e-8)}}));
  monte_carlo->attempt(1e2);
  monte_carlo->add(MakeNumParticles({{"steps_per_write", str(steps_per)},
                                    {"file_name", "tmp/pljavbnum.txt"}}));
  monte_carlo->add(MakeEnergy({{"steps_per_write", str(steps_per)},
                              {"file_name", "tmp/pljavbe.txt"}}));
  monte_carlo->attempt(1e2);
}

}  // namespace feasst
