#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "utils/include/progress_report.h"
#include "math/include/random_mt19937.h"
#include "system/include/visit_model.h"
#include "system/include/visit_model_inner.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "system/include/lennard_jones.h"
#include "models/include/square_well.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/remove.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/tune.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/energy.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
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

///// This test was decommissioned on 12/3/2024 - energy errors upon deserialization or CheckEnergy updates > 1
///// Without single particle translations, rigid cluster moves should reject
///// cluster coalescence and breakup to satisfy detailed balance.
//TEST(MonteCarlo, cluster_LONG) {
//  //for (auto single_particle_translate : {false}) {
//  //for (auto single_particle_translate : {true}) {
//  for (auto single_particle_translate : {true, false}) {
//    MonteCarlo mc;
//    //mc.set(MakeRandomMT19937({{"seed", "1613161559"}}));
//    mc.add(MakeConfiguration({{"cubic_side_length", "8"},
//      {"particle_type0", "../particle/lj.txt"},
//      {"add_particles_of_type0", "3"}}));
//    mc.get_system()->get_configuration()->update_positions({{0, 0, 0}, {2, 0, 0}, {4, 0, 0}});
//    mc.add(MakePotential({{"Model", "LennardJones"}, {"EnergyMap", "EnergyMapNeighbor"}}));
//    //mc.add(MakePotential({{"Model", "LennardJones"}, {"VisitModel", "VisitModel"}, {"VisitModelInner", "VisitModelInner"}, {"EnergyMap", "EnergyMapNeighbor"}}));
//    //mc.add(MakePotential(MakeLennardJones(),
//    //  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighbor()))));
//      //MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
//    mc.set(MakeThermoParams({{"beta", "40"}, {"chemical_potential", "1."}}));
//    mc.set(MakeMetropolis());
//    if (single_particle_translate) mc.add(MakeTrialTranslate());
//    mc.add(MakeNeighborCriteria({{"energy_maximum", "-0.5"}}));
//    auto scluster = MakeSelectCluster({{"neighbor_index", "0"}});
//    scluster->select_cluster(0, mc.system());
//    const int cluster_size = scluster->mobile().num_particles();
//    mc.add(MakeTrialRigidCluster({
//      {"particle_type", "0"},
//      {"neighbor_index", "0"},
//      {"rotate_param", "50"},
//      {"translate_param", "1"}}));
//    const int trials_per = 1e2;
////    mc.add(MakeLogAndMovie({{"trials_per_write", str(trials_per)}, {"output_file", "tmp/lj"}}));
//    mc.add(MakeCheckEnergy({{"trials_per_update", "1e1"}}));
//    mc.add(MakeTune({{"trials_per_write", str(trials_per)}, {"output_file", "tmp/tune.txt"}}));
//    // conduct the trials
//    const VisitModelInner& inner = mc.system().potential(0).visit_model().inner();
//    for (int trial = 0; trial < 1e5; ++trial) {
//      //INFO("trial " << trial);
//      mc.attempt(1);
//      EXPECT_NEAR(inner.energy_map().total_energy(),
//                  mc.system().stored_energy(),
//                  NEAR_ZERO);
//    }
//    scluster->select_cluster(0, mc.system());
//    if (single_particle_translate) {
//      EXPECT_EQ(scluster->mobile().num_particles(),
//                mc.configuration().num_particles());
//    } else {
//      EXPECT_EQ(scluster->mobile().num_particles(), cluster_size);
//    }
//    auto mc2 = test_serialize_unique(mc);
////    mc2->attempt(1e1);
////    // ensure TrialFactory still tunes
////    if (single_particle_translate) {
////      int rigid_index = 1;
////      EXPECT_NE(mc.trial(rigid_index).stage(0).perturb().tunable().value(), 1.);
////      EXPECT_NE(mc.trial(rigid_index+1).stage(0).perturb().tunable().value(), 50.);
////    }
//  }
//}

  //  // test rotation through PBCs
  //  mc.get_system()->get_configuration()->update_positions({
  //    {4.0, 0.769800358919501000, 0},
  //    { 4.666666666666666666, -0.384900179459751000, 0},
  //    { 3.333333333333333334, -0.384900179459751000, 0}});
//    mc.get_system()->energy();
//    mc.add(MakeMetropolis({{"beta", "40"}, {"chemical_potential", "1."}}));

TEST(MonteCarlo, GCMCmap) {
  //for (std::string mapstr : {"neighbor"}) {
  for (std::string mapstr : {"all", "neighbor"}) {
    INFO(mapstr);
    MonteCarlo mc;
    mc.set(MakeRandomMT19937({{"seed", "123"}}));
    mc.add(MakeConfiguration({{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}));
    mc.add(MakePotential(MakeLennardJones()));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    mc.set(MakeMetropolis());
    mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
    const std::string trials_per = str(1e4);
//    mc.add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/lj"}}));
    mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}, {"decimal_places", "6"}}));
    mc.add(MakeTune());
    std::shared_ptr<EnergyMap> map;
    if (mapstr == "all") {
      map = MakeEnergyMapAll();
    } else if (mapstr == "neighbor") {
      map = MakeEnergyMapNeighbor();
    } else {
      FATAL("unrecognized mapstr");
    }
    mc.set(0, MakePotential(MakeLennardJones(),
      MakeVisitModel(MakeVisitModelInner(map))));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "-4"}}));
    mc.add(MakeTrialTransfer({{"particle_type", "0"}}));
    mc.add(MakeNumParticles({{"trials_per_write", trials_per},
                             {"output_file", "tmp/ljnum.txt"}}));
    for (int i = 0; i < 1e4; ++i) {
      mc.attempt(1);
      const double en = mc.criteria().current_energy();
      const double en_map = mc.system().potential(0).visit_model().inner().energy_map().total_energy();
      if (std::abs(en - en_map) > 1e-8) {
        FATAL(MAX_PRECISION << "not the same: " << en << " " << en_map);
      }
    }
  }
}

std::unique_ptr<MonteCarlo> mc_avb_test(
    const bool avb = true,
    const int min_particles = 1,
    const bool avb2 = false,
    const bool avb4 = false) {
  auto monte_carlo = std::make_unique<MonteCarlo>();
  // monte_carlo->set(MakeRandomMT19937({{"seed", "default"}}));
  monte_carlo->add(MakeConfiguration({{"cubic_side_length", "6"},
                                     {"particle_type", "../particle/lj.txt"}}));
  if (avb) {
    monte_carlo->add(MakePotential({{"Model", "LennardJones"}, {"EnergyMap", "EnergyMapAll"}}));
    //monte_carlo->add(MakePotential(MakeLennardJones(),
    //  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
      //MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighbor()))));
  } else {
    monte_carlo->add(MakePotential(MakeLennardJones()));
  }
  monte_carlo->set(MakeThermoParams({{"beta", "0.00001"}, {"chemical_potential", "50."}}));
  monte_carlo->set(MakeMetropolis());
  monte_carlo->add(MakeTrialAdd({{"particle_type", "0"}}));
  monte_carlo->run(MakeRun({{"until_num_particles", str(min_particles)}}));
  monte_carlo->run(MakeRemove({{"name", "TrialAdd"}}));
  monte_carlo->set(MakeThermoParams({{"beta", "0.2"}, {"chemical_potential", "-20."}}));
  monte_carlo->set(MakeMetropolis(
    MakeConstrainNumParticles({{"minimum", str(min_particles)}})));
  if (avb) {
    monte_carlo->add(MakeNeighborCriteria({{"maximum_distance", "3"},
                                          {"minimum_distance", "1"}}));
    if (avb2) {
      monte_carlo->add(MakeTrialAVB2({{"neighbor_index", "0"},
                                     {"particle_type", "0"}}));
    } else if (avb4) {
      monte_carlo->add(MakeTrialAVB4({{"neighbor_index", "0"},
                                     {"particle_type", "0"}}));
    } else {
      monte_carlo->add(MakeTrialTransferAVB(
        {{"particle_type", "0"}, {"neighbor_index", "0"}}));
    }
  } else {
    if (avb2 || avb4) {
      monte_carlo->add(MakeTrialTranslate({{"tunable_param", "4"}}));
    } else {
      monte_carlo->add(MakeTrialTransfer({{"particle_type", "0"}}));
    }
  }
  const int trials_per = 1e4;
  monte_carlo->add(MakeMovie({{"trials_per_write", str(trials_per)},
                             {"output_file", "tmp/ljavb.xyz"},
                             {"clear_file", "true"}}));
  monte_carlo->add(MakeLog({{"trials_per_write", str(trials_per)},
                           {"output_file", "tmp/ljavb_log.txt"},
                           {"clear_file", "true"}}));
  monte_carlo->add(MakeCheckEnergy({{"trials_per_update", str(trials_per)},
                                   {"tolerance", str(1e-8)}}));
  monte_carlo->attempt(1e6);
  monte_carlo->add(MakeNumParticles({{"trials_per_write", str(trials_per)},
                                    {"output_file", "tmp/ljavbnum.txt"}}));
  monte_carlo->add(MakeEnergy({{"trials_per_write", str(trials_per)},
                              {"output_file", "tmp/ljavbe.txt"}}));
  monte_carlo->attempt(1e6);
  return test_serialize_unique(*monte_carlo);
}

const double z_factor = 20.;

TEST(MonteCarlo, GCMC_AVB_LONG) {
  auto mc_avb = mc_avb_test(true);
  auto mc_noavb = mc_avb_test(false);
  INFO(mc_avb->analyze(2).accumulator().str())
  INFO(mc_avb->analyze(3).accumulator().str())
  INFO(mc_noavb->analyze(2).accumulator().str())
  INFO(mc_noavb->analyze(3).accumulator().str())
  EXPECT_TRUE(mc_avb->analyze(2).accumulator().is_equivalent(
            mc_noavb->analyze(2).accumulator(), z_factor, 0, true));
//  EXPECT_TRUE(mc_avb->analyze(3).accumulator().is_equivalent(
//            mc_noavb->analyze(3).accumulator(), 3, 0, true));
}

TEST(MonteCarlo, MC_AVB2_AVB4_LONG) {
  const int num_particles = 10;
  auto mc_noavb = mc_avb_test(false, num_particles, true);
  EXPECT_NEAR(10, mc_noavb->analyze(2).accumulator().average(), NEAR_ZERO);
  INFO(mc_noavb->analyze(3).accumulator().str())

  auto mc_avb2 = mc_avb_test(true, num_particles, true);
  INFO(mc_avb2->analyze(3).accumulator().str())
  EXPECT_NEAR(10, mc_avb2->analyze(2).accumulator().average(), NEAR_ZERO);
  EXPECT_TRUE(mc_avb2->analyze(3).accumulator().is_equivalent(
             mc_noavb->analyze(3).accumulator(), z_factor, 0, true));

  if (true) {
    auto mc_avb4 = mc_avb_test(true, num_particles, false, true);
    EXPECT_NEAR(10, mc_avb4->analyze(2).accumulator().average(), NEAR_ZERO);
    INFO(mc_avb4->analyze(3).accumulator().str())
    EXPECT_TRUE(mc_avb4->analyze(3).accumulator().is_equivalent(
               mc_noavb->analyze(3).accumulator(), z_factor, 0, true));
  }
}

TEST(MonteCarlo, binaryavb) {
  auto mc = std::make_unique<MonteCarlo>();
  mc->set(MakeRandomMT19937({{"seed", "12345"}}));
  mc->add(MakeConfiguration({
    {"xyz_file", "../plugin/cluster/test/data/ab.xyz"},
    {"particle_type", "A:../plugin/cluster/test/data/a.txt,B:../plugin/cluster/test/data/b.txt"}}));
  mc->add(MakeNeighborCriteria({{"site_type0", "A,A,B"}, {"site_type1", "A,B,B"}, {"maximum_distance", "2.075,1.575,1.075"}, {"minimum_distance", "2,1.5,1"}, {"energy_maximum", "1e9,1e9,1e9"}}));
  mc->add(MakeNeighborCriteria({{"site_type0", "A"}, {"site_type1", "A"}, {"maximum_distance", "2.075"}, {"minimum_distance", "2"}}));
  mc->add(MakeNeighborCriteria({{"site_type0", "A"}, {"site_type1", "B"}, {"maximum_distance", "1.575"}, {"minimum_distance", "1.5"}}));
  mc->add(MakeNeighborCriteria({{"site_type0", "B"}, {"site_type1", "B"}, {"maximum_distance", "1.075"}, {"minimum_distance", "1"}}));
  mc->add(MakePotential({{"Model", "LennardJones"}, {"EnergyMap", "EnergyMapNeighborCriteria"}, {"neighbor_index", "0"}}));
  mc->set(MakeThermoParams({{"beta", "0.2"}, {"chemical_potential", "-20."}}));
  mc->set(MakeMetropolis());
  //mc->add(MakeTrialAVB2({{"particle_type", "A"}, {"site", "A0"}, {"target_particle_type", "A"}, {"target_site", "A0"}, {"neighbor_index", "1"}}));
  //mc->add(MakeTrialAVB2({{"particle_type", "A"}, {"site", "A0"}, {"target_particle_type", "B"}, {"target_site", "B0"}, {"neighbor_index", "2"}}));
  mc->add(MakeTrialAVB2({{"particle_type", "B"}, {"site", "B0"}, {"target_particle_type", "A"}, {"target_site", "A0"}, {"neighbor_index", "2"}}));
  //mc->add(MakeTrialAVB2({{"particle_type", "B"}, {"site", "B0"}, {"target_particle_type", "B"}, {"target_site", "B0"}, {"neighbor_index", "3"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", "1"}, {"decimal_places", "6"}}));
  //mc->attempt(1);
  auto mc2 = test_serialize_unique(*mc);
  EXPECT_EQ(mc2->configuration().particle(0).type(), 0);
  EXPECT_EQ(mc2->configuration().particle(1).type(), 1);
  mc2->attempt(2);
  EXPECT_EQ(mc2->trials().num_success(), 1);
  EXPECT_NEAR(mc2->trial(0).accept().ln_metropolis_prob(), -6.981754116050470000, 1e-13);
  //INFO("U_new:" << MAX_PRECISION << mc2->trials().trial(0).accept().energy_new());
}

std::unique_ptr<MonteCarlo> binary_avb_test(const bool avb) {
  auto mc = std::make_unique<MonteCarlo>();
//  mc->set(MakeRandomMT19937({{"seed", "12345"}}));
  mc->add(MakeConfiguration({
    {"particle_type", "A:../plugin/cluster/test/data/a.txt,B:../plugin/cluster/test/data/b.txt"},
    {"xyz_file", "../plugin/cluster/test/data/ab5.xyz"}}));
  auto nc0 = MakeNeighborCriteria({{"site_type0", "A,A,B"}, {"site_type1", "A,B,B"}, {"maximum_distance", "5.075,3.575,2.075"}, {"minimum_distance", "2,1.5,1"}, {"energy_maximum", "1e9,1e9,1e9"}});
  mc->add(nc0);
  auto nc1 = MakeNeighborCriteria({{"site_type0", "A"}, {"site_type1", "A"}, {"maximum_distance", "5.075"}, {"minimum_distance", "2"}});
  mc->add(nc1);
  auto nc2 = MakeNeighborCriteria({{"site_type0", "A"}, {"site_type1", "B"}, {"maximum_distance", "3.575"}, {"minimum_distance", "1.5"}});
  mc->add(nc2);
  auto nc3 = MakeNeighborCriteria({{"site_type0", "B"}, {"site_type1", "B"}, {"maximum_distance", "2.075"}, {"minimum_distance", "1"}});
  mc->add(nc3);
  //EXPECT_EQ(mc->system().neighbor_criteria(0, 0).site_type0(), 0);
  //EXPECT_EQ(mc->system().neighbor_criteria(0, 0).site_type1(), 0);
  EXPECT_EQ(mc->system().neighbor_criteria(1, 0).site_type0(), 0);
  EXPECT_EQ(mc->system().neighbor_criteria(1, 0).site_type1(), 0);
  EXPECT_EQ(mc->system().neighbor_criteria(2, 0).site_type0(), 0);
  EXPECT_EQ(mc->system().neighbor_criteria(2, 0).site_type1(), 1);
  EXPECT_EQ(mc->system().neighbor_criteria(3, 0).site_type0(), 1);
  EXPECT_EQ(mc->system().neighbor_criteria(3, 0).site_type1(), 1);
  mc->add(MakePotential({{"Model", "LennardJones"}, {"EnergyMap", "EnergyMapNeighborCriteria"}, {"neighbor_index", "0"}}));
  mc->set(MakeThermoParams({{"beta", "0.2"}, {"chemical_potential", "-20."}}));
  mc->set(MakeMetropolis());
  if (avb) {
    mc->add(MakeTrialAVB2({{"particle_type", "A"}, {"site", "A0"}, {"target_particle_type", "A"}, {"target_site", "A0"}, {"neighbor_index", "1"}}));
    mc->add(MakeTrialAVB2({{"particle_type", "A"}, {"site", "A0"}, {"target_particle_type", "B"}, {"target_site", "B0"}, {"neighbor_index", "2"}}));
    mc->add(MakeTrialAVB2({{"particle_type", "B"}, {"site", "B0"}, {"target_particle_type", "A"}, {"target_site", "A0"}, {"neighbor_index", "2"}}));
    mc->add(MakeTrialAVB2({{"particle_type", "B"}, {"site", "B0"}, {"target_particle_type", "B"}, {"target_site", "B0"}, {"neighbor_index", "3"}}));
    mc->add(MakeTrialAVB4({{"particle_type", "A"}, {"site", "A0"}, {"target_particle_type", "A"}, {"target_site", "A0"}, {"neighbor_index", "1"}}));
    mc->add(MakeTrialAVB4({{"particle_type", "A"}, {"site", "A0"}, {"target_particle_type", "B"}, {"target_site", "B0"}, {"neighbor_index", "2"}}));
    mc->add(MakeTrialAVB4({{"particle_type", "B"}, {"site", "B0"}, {"target_particle_type", "A"}, {"target_site", "A0"}, {"neighbor_index", "2"}}));
    mc->add(MakeTrialAVB4({{"particle_type", "B"}, {"site", "B0"}, {"target_particle_type", "B"}, {"target_site", "B0"}, {"neighbor_index", "3"}}));
  } else {
    mc->add(MakeTrialTranslate());
  }
  mc->add(MakeCheckEnergy({{"trials_per_update", "1e4"}, {"decimal_places", "6"}}));
  mc->add(MakeEnergy({{"trials_per_write", "1e4"}, {"output_file", "tmp/bavben.csv"}}));
  auto mc2 = test_serialize_unique(*mc);
  mc2->attempt(1e6);
  return mc2;
}

TEST(MonteCarlo, binaryavb_LONG) {
  auto mc_avb = binary_avb_test(true);
  DEBUG(mc_avb->analyze(0).accumulator().str());
  auto mc_noavb = binary_avb_test(false);
  DEBUG(mc_noavb->analyze(0).accumulator().str());
  EXPECT_TRUE(mc_avb->analyze(0).accumulator().is_equivalent(
            mc_noavb->analyze(0).accumulator(), z_factor, 0, true));
}

}  // namespace feasst
