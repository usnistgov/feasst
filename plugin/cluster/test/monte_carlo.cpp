#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "system/include/utils.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/energy.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_avb4.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/trial_transfer_avb.h"
//#include "cluster/include/trial_add_avb.h"
//#include "cluster/include/trial_remove_avb.h"

namespace feasst {

/// Without single particle translations, rigid cluster moves should reject
/// cluster coalescence and breakup to satisfy detailed balance.
TEST(MonteCarlo, cluster) {
  for (auto single_particle_translate : {true, false}) {
    MonteCarlo mc;
    // mc.set(MakeRandomMT19937({{"seed", "default"}}));
    // mc.set(MakeRandomMT19937({{"seed", "1580855528"}}));
    mc.add(Configuration(MakeDomain({{"cubic_box_length", "8"}}),
                                  {{"particle_type", "../forcefield/data.lj"}}));
    mc.add(Potential(MakeLennardJones(),
      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
    mc.add(MakeMetropolis({{"beta", "40"}, {"chemical_potential", "1."}}));
    if (single_particle_translate) mc.add(MakeTrialTranslate());
    SeekNumParticles(3).with_trial_add().run(&mc);
    auto neighbor_criteria = MakeNeighborCriteria({{"energy_maximum", "-0.5"}});
    SelectCluster scluster(neighbor_criteria);
    scluster.select_cluster(0, mc.system());
    const int cluster_size = scluster.mobile().num_particles();
    mc.add(MakeTrialRigidCluster(
      neighbor_criteria,
      { {"rotate_param", "50"},
        {"translate_param", "1"}}));
    const int steps_per = 1e0;
    mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/lj"}}));
    mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
    // conduct the trials
    const VisitModelInner& inner = mc.system().potential(0).visit_model().inner();
    for (int trial = 0; trial < 1e3; ++trial) {
      //INFO("trial " << trial);
      mc.attempt(1);
      EXPECT_NEAR(inner.energy_map().total_energy(),
                  mc.system().stored_energy(),
                  NEAR_ZERO);
    }
    scluster.select_cluster(0, mc.system());
    if (single_particle_translate) {
      EXPECT_EQ(scluster.mobile().num_particles(),
                mc.configuration().num_particles());
    } else {
      EXPECT_EQ(scluster.mobile().num_particles(), cluster_size);
    }
  }
}

  //  // test rotation through PBCs
  //  mc.get_system()->get_configuration()->update_positions({
  //    {4.0, 0.769800358919501000, 0},
  //    { 4.666666666666666666, -0.384900179459751000, 0},
  //    { 3.333333333333333334, -0.384900179459751000, 0}});
//    mc.get_system()->energy();
//    mc.add(MakeMetropolis({{"beta", "40"}, {"chemical_potential", "1."}}));

TEST(MonteCarlo, GCMCmap) {
  MonteCarlo mc;
  mc.set(lennard_jones({{"lrc", "false"}}));
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}}));
  mc.set(0, Potential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-6"}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}}));
  mc.add(MakeNumParticles({{"steps_per_write", "1000"},
                           {"file_name", "tmp/ljnum.txt"}}));
  for (int i = 0; i < 1e4; ++i) {
    mc.attempt(1);
    const double en = mc.criteria().current_energy();
    const double en_map = mc.system().potential(0).visit_model().inner().energy_map().total_energy();
    if (std::abs(en - en_map) > 1e-8) {
      INFO(MAX_PRECISION << "not the same: " << en << " " << en_map);
    }
  }
}

MonteCarlo mc_avb_test(
    const bool avb = true,
    const int min_particles = 1,
    const bool avb2 = false,
    const bool avb4 = false) {
  MonteCarlo monte_carlo;
  // monte_carlo.set(MakeRandomMT19937({{"seed", "default"}}));
  monte_carlo.add(Configuration(MakeDomain({{"cubic_box_length", "6"}}),
                                {{"particle_type", "../forcefield/data.lj"}}));
  if (avb) {
    monte_carlo.add(Potential(MakeLennardJones(),
      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
      //MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighbor()))));
  } else {
    monte_carlo.add(Potential(MakeLennardJones()));
  }
  monte_carlo.add(MakeMetropolis({{"beta", "0.00001"}, {"chemical_potential", "50."}}));
  SeekNumParticles(min_particles).with_trial_add().run(&monte_carlo);
  monte_carlo.add(MakeMetropolis(
    MakeConstrainNumParticles({{"minimum", str(min_particles)}}),
    {{"beta", "0.2"}, {"chemical_potential", "-20."}}));
  if (avb) {
    auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "3"},
                                                   {"minimum_distance", "1"}});
    if (avb2) {
      monte_carlo.add(MakeTrialAVB2(neighbor_criteria));
    } else if (avb4) {
      monte_carlo.add(MakeTrialAVB4(neighbor_criteria));
    } else {
      monte_carlo.add(MakeTrialTransferAVB(neighbor_criteria,
        {{"particle_type", "0"}}));
    }
  } else {
    if (avb2 || avb4) {
      monte_carlo.add(MakeTrialTranslate({{"tunable_param", "4"}}));
    } else {
      monte_carlo.add(MakeTrialTransfer({{"particle_type", "0"}}));
    }
  }
  const int steps_per = 1e4;
  monte_carlo.add(MakeMovie({{"steps_per", str(steps_per)},
                             {"file_name", "tmp/ljavb.xyz"},
                             {"clear_file", "true"}}));
  monte_carlo.add(MakeLog({{"steps_per", str(steps_per)},
                           {"file_name", "tmp/ljavb_log.txt"},
                           {"clear_file", "true"}}));
  monte_carlo.add(MakeCheckEnergy({{"steps_per", str(steps_per)},
                                   {"tolerance", str(1e-8)}}));
  monte_carlo.attempt(1e6);
  monte_carlo.add(MakeNumParticles({{"steps_per_write", str(steps_per)},
                                    {"file_name", "tmp/ljavbnum.txt"}}));
  monte_carlo.add(MakeEnergy({{"steps_per_write", str(steps_per)},
                              {"file_name", "tmp/ljavbe.txt"}}));
  monte_carlo.attempt(1e6);
  return monte_carlo;
}

const double z_factor = 10.;

TEST(MonteCarlo, GCMC_AVB_LONG) {
  MonteCarlo mc_avb = mc_avb_test(true);
  MonteCarlo mc_noavb = mc_avb_test(false);
  INFO(mc_avb.analyze(2).accumulator().str())
  INFO(mc_avb.analyze(3).accumulator().str())
  INFO(mc_noavb.analyze(2).accumulator().str())
  INFO(mc_noavb.analyze(3).accumulator().str())
  EXPECT_TRUE(mc_avb.analyze(2).accumulator().is_equivalent(
            mc_noavb.analyze(2).accumulator(), z_factor, true));
//  EXPECT_TRUE(mc_avb.analyze(3).accumulator().is_equivalent(
//            mc_noavb.analyze(3).accumulator(), 3, true));
}

TEST(MonteCarlo, MC_AVB2_AVB4_LONG) {
  const int num_particles = 10;
  MonteCarlo mc_noavb = mc_avb_test(false, num_particles, true);
  EXPECT_NEAR(10, mc_noavb.analyze(2).accumulator().average(), NEAR_ZERO);
  INFO(mc_noavb.analyze(3).accumulator().str())

  MonteCarlo mc_avb2 = mc_avb_test(true, num_particles, true);
  INFO(mc_avb2.analyze(3).accumulator().str())
  EXPECT_NEAR(10, mc_avb2.analyze(2).accumulator().average(), NEAR_ZERO);
  EXPECT_TRUE(mc_avb2.analyze(3).accumulator().is_equivalent(
             mc_noavb.analyze(3).accumulator(), z_factor, true));

  if (true) {
    MonteCarlo mc_avb4 = mc_avb_test(true, num_particles, false, true);
    EXPECT_NEAR(10, mc_avb4.analyze(2).accumulator().average(), NEAR_ZERO);
    INFO(mc_avb4.analyze(3).accumulator().str())
    EXPECT_TRUE(mc_avb4.analyze(3).accumulator().is_equivalent(
               mc_noavb.analyze(3).accumulator(), z_factor, true));
  }
}

}  // namespace feasst
