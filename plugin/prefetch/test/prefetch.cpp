#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
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
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_avb4.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/trial_transfer_avb.h"

namespace feasst {

void run_prefetch(const int trials, const int trials_per) {
  auto mc = MakePrefetch();
//  mc->set(MakeRandomMT19937({{"seed", "1592943710"}}));
  mc->set(MakeRandomMT19937({{"seed", "1596650884"}}));
  mc->add(MakeConfiguration({{"cubic_box_length", "8"},
                             {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc->add(MakePotential(MakeLennardJones()));
  mc->add(MakePotential(MakeLongRangeCorrections()));
  mc->set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}, {"num_steps", "1"}}));
  mc->add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/lj"}}));
  mc->add(MakeCheckEnergyAndTune({{"trials_per", str(trials_per)}}));
  mc->activate_prefetch(false);
  mc->add(MakeTrialAdd({{"particle_type", "0"}}));
  mc->run(MakeRun({{"until_num_particles", "50"}}));
  mc->run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  // activate prefetch after initial configuration
  mc->activate_prefetch(true);
  mc->attempt(trials);
  EXPECT_EQ(mc->analyze(0).trials_since_write(),
            mc->modify(0).trials_since_update());
}

TEST(Prefetch, NVT_benchmark) {
  run_prefetch(1e3, 1e1);
}

TEST(Prefetch, NVT_benchmark_LONG) {
  run_prefetch(1e6, 1e3); // 5.4s on 4 cores of i7-4770K @ 3.5GHz
}

void prefetch(System system, const int sync = 0) {
  auto mc = MakePrefetch({{"trials_per_check", "1"}, {"synchronize", str(sync)}});
  //mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->set(MakeRandomMT19937({{"seed", "time"}}));
  mc->set(system);
  mc->set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeLogAndMovie({{"trials_per", str(1e1)}, {"file_name", "tmp/lj"}}));
  mc->add(MakeCheckEnergyAndTune({{"trials_per", str(1e1)}}));
  //mc_lj(mc.get(), 8, "../forcefield/lj.fstprt", 1e1, true, false);
  // mc->set(MakeRandomMT19937({{"seed", "default"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578665877"}}));
  // mc->set(MakeRandomMT19937({{"seed", "1578667496"}}));
  //mc->set(MakeRandomMT19937({{"seed", "1804289383"}}));
  mc->add(MakeTrialAdd({{"particle_type", "0"}}));
  mc->add(MakeTrialRemove({{"particle_type", "0"}}));
  // mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-2"}}));
  mc->set(MakeThermoParams({{"beta", str(1./1.5)},
     {"chemical_potential", "-2.352321"}}));
  mc->set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}})),
    MakeTransitionMatrix({{"min_sweeps", "10"}})));
  mc->add(MakeCriteriaUpdater({{"trials_per", str(1e1)}}));
  mc->activate_prefetch(true);
  mc->attempt(10);
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
  EXPECT_EQ(1, mc2.trials_per_check());
}

TEST(Prefetch, MUVT) {
  System sys;
  sys.add(MakeConfiguration({{"cubic_box_length", "8"},
                             {"particle_type0", "../forcefield/lj.fstprt"}}));
  sys.add(MakePotential(MakeLennardJones()));
  sys.add(MakePotential(MakeLongRangeCorrections()));
  prefetch(sys);
}

TEST(Prefetch, MUVT_spce) {
  prefetch(spce({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}, {"table_size", str(1e3)}}), 1);
}

TEST(Prefetch, NVT_spce) {
  auto mc = MakePrefetch({{"synchronize", "true"}});
  //auto mc = MakePrefetch({{"synchronize", "false"}});
  // mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->set(spce({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}, {"table_size", str(1e3)}}));
  const int trials_per = 1e2;
  mc->set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/lj"}}));
  //mc->add(MakeCheckProperties({{"trials_per", "1"}}));
  mc->add(MakeCheckProperties({{"trials_per", str(trials_per)}}));
  mc->add(MakeCheckEnergyAndTune({{"trials_per", str(trials_per)}}));
  // mc->set(MakeRandomMT19937({{"seed", "default"}}));
  mc->activate_prefetch(false);
  mc->add(MakeTrialAdd({{"particle_type", "0"}}));
  mc->run(MakeRun({{"until_num_particles", "50"}}));
  mc->run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  // activate prefetch after initial configuration
  mc->activate_prefetch(true);
  // mc->attempt(1e6);  // ~3.5 seconds (now 4.1)
  mc->attempt(1e2);
  EXPECT_EQ(mc->analyze(0).trials_since_write(),
            mc->modify(1).trials_since_update());
}

TEST(Prefetch, AVB) {
  const int trials_per = 1e2;
  auto monte_carlo = MakePrefetch({{"synchronize", "true"},
                                   {"trials_per_check", str(trials_per)}});
  monte_carlo->add(MakeConfiguration({{"cubic_box_length", "6"},
                                      {"particle_type", "../forcefield/lj.fstprt"}}));
  monte_carlo->add(MakePotential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  //monte_carlo->add(MakePotential(MakeLennardJones()));
  monte_carlo->set(MakeThermoParams({{"beta", "0.00001"}, {"chemical_potential", "50."}}));
  monte_carlo->set(MakeMetropolis());
  monte_carlo->activate_prefetch(false);
  monte_carlo->add(MakeTrialAdd({{"particle_type", "0"}}));
  monte_carlo->run(MakeRun({{"until_num_particles", "50"}}));
  monte_carlo->run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  monte_carlo->activate_prefetch(true);
  monte_carlo->set(MakeThermoParams({{"beta", "0.2"}, {"chemical_potential", "-20."}}));
  monte_carlo->add(MakeNeighborCriteria({{"maximum_distance", "3"},
                                         {"minimum_distance", "1"}}));
  // Something wrong with adding TrialFactories when using prefetch...
  // I think they don't end up selecting the right trial type
//  monte_carlo->add(MakeTrialAVB2({{"neighbor_index", "0"}}));
  monte_carlo->add(MakeTrialAVB2Half({{"particle_type", "0"}, {"neighbor_index", "0"}, {"out_to_in", "true"}}));
  monte_carlo->add(MakeTrialAVB2Half({{"particle_type", "0"}, {"neighbor_index", "0"}, {"out_to_in", "false"}}));
  monte_carlo->add(MakeTrialAVB4({{"particle_type", "0"}, {"neighbor_index", "0"}}));
//  monte_carlo->add(MakeTrialTransferAVB({{"neighbor_index", "0"},
//    {"particle_type", "0"}}));
  monte_carlo->add(MakeTrialTranslate({{"tunable_param", "4"}}));
//  monte_carlo->add(MakeTrialTransfer({{"particle_type", "0"}}));
  monte_carlo->add(MakeLogAndMovie({{"trials_per", str(trials_per)},
                             {"file_name", "tmp/pljavb.xyz"},
                             {"clear_file", "true"}}));
  monte_carlo->add(MakeCheckEnergy({{"trials_per", str(trials_per)},
                                   {"tolerance", str(1e-8)}}));
  monte_carlo->attempt(1e2);
  monte_carlo->add(MakeNumParticles({{"trials_per_write", str(trials_per)},
                                    {"file_name", "tmp/pljavbnum.txt"}}));
  monte_carlo->add(MakeEnergy({{"trials_per_write", str(trials_per)},
                              {"file_name", "tmp/pljavbe.txt"}}));
  monte_carlo->attempt(1e2);
}

}  // namespace feasst
