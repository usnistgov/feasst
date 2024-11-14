#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "configuration/include/physical_constants.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/tabulate_two_rigid_body_3D.h"

namespace feasst {

inline void run_hs(const int num_orientations_per_pi, const int proc, const int num_proc) {
  auto mc_hs = MakeMonteCarlo({{
    {"Configuration", {
      {"cubic_side_length", "2e2"},
      {"particle_type0", "../plugin/aniso/test/data/fc.fstprt"},
      {"particle_type1", "../plugin/aniso/test/data/fc.fstprt"},
      {"add_particles_of_type0", "1"},
      {"add_particles_of_type1", "1"},
      {"group0", "fixed"},
      {"fixed_particle_type", "0"},
      {"group1", "mobile"},
      {"mobile_particle_type", "1"},
      {"set_cutoff_min_to_sigma", "true"},
    }},
    {"Potential", {{"Model", "HardSphere"}, {"VisitModel", "VisitModelCell"}, {"min_length", "max_sigma"}, {"energy_cutoff", "1e100"}}},
  }});

//  INFO("fraction unique " << table_ior->rotator().fraction_unique());
  auto table_hs = std::make_shared<TabulateTwoRigidBody3D>(argtype({
    {"proc", str(proc)}, {"num_proc", str(num_proc)},
    {"input_orientation_file", "tmp/unique_ior.txt"},
    {"num_z", "-1"},
    {"output_table_file", "tmp/contact_table"+str(proc)+".txt"},
  }));
  EXPECT_NEAR(207.84639874816824, table_hs->max_cubic_side_length(0, mc_hs->configuration()), NEAR_ZERO);
  table_hs->run(mc_hs.get());
  EXPECT_NEAR(1./9., table_hs->rotator().fraction_unique(), NEAR_ZERO);
  EXPECT_NEAR(8978990.5049453303, mc_hs->configuration().domain().volume(), 1e-6);
  if (num_proc == 1) {
    EXPECT_EQ(72, table_hs->rotator().num_orientations());
  } else if (num_proc == 2) {
    EXPECT_EQ(36, table_hs->rotator().num_orientations());
  }
  int ior = 0;
  for (int iorall = 0; iorall < 3; ++iorall) {
    if (table_hs->rotator().ior_in_proc(iorall)) {
      double expect;
      if (iorall == 0 || iorall == 2) expect = 74.09857177734375;
      //if (iorall == 0 || iorall == 2) expect = 74.098304487805038;
      if (iorall == 1) expect = 71.761772155761719;
      //if (iorall == 1) expect = 71.761553822326817;
      EXPECT_NEAR(expect, table_hs->rotator().contact_[ior], 5e-5);
      ++ior;
    }
  }
  // generate single contact file
  if (proc == num_proc - 1) {
    table_hs->combine_table({{"prefix", "tmp/contact_table"},
                             {"suffix", ".txt"}});
  }
}

inline void run(const int num_orientations_per_pi, const int proc, const int num_proc) {
  const double temperature = 298.15;
  const double ionic_strength = 0.15;
  const double smoothing_distance = 2.;
  const double temp_cel = temperature - 273.15;
  const double dielectric_water = 87.74 - 0.40008*temp_cel + 9.398e-4*std::pow(temp_cel, 2) - 1.4e-6*std::pow(temp_cel, 3);
  EXPECT_NEAR(dielectric_water, 78.3035, 1e-4);
  const double eps_0 = CODATA2018().permitivity_vacuum();
  const double elem_q = CODATA2018().elementary_charge();
  const double na = CODATA2018().avogadro_constant();
  const double kb = CODATA2018().boltzmann_constant();
  const double kappa = std::sqrt(2*std::pow(elem_q, 2)*ionic_strength*(1e3)*na/(dielectric_water*eps_0*kb*temperature*1e20));
  EXPECT_NEAR(kappa, 0.1274742518905759, NEAR_ZERO);
  const double cutoff = 5./kappa;
  EXPECT_NEAR(cutoff, 39.22360732339899, NEAR_ZERO);

  auto mc = MakeMonteCarlo({{
    {"Configuration", {
      {"cubic_side_length", "2e2"},
      {"particle_type0", "../plugin/aniso/test/data/fc.fstprt"},
      {"particle_type1", "../plugin/aniso/test/data/fc.fstprt"},
      {"add_particles_of_type0", "1"},
      {"add_particles_of_type1", "1"},
      {"group0", "fixed"},
      {"fixed_particle_type", "0"},
      {"group1", "mobile"},
      {"mobile_particle_type", "1"},
      {"cutoff", str(cutoff)},
    }},
    {"Potential", {{"Model", "ModelTwoBodyFactory"},
      {"model0", "HardSphere"},
      {"model1", "LennardJones"},
      {"model2", "DebyeHuckel"}, {"kappa", str(kappa)}, {"dielectric", str(dielectric_water)}, {"smoothing_distance", str(smoothing_distance)},
      {"VisitModel", "VisitModelCell"}, {"min_length", "max_cutoff"}, {"energy_cutoff", "1e100"}}},
  }});
  auto table = std::make_shared<TabulateTwoRigidBody3D>(argtype({
    {"proc", str(proc)}, {"num_proc", str(num_proc)},
    {"input_orientation_file", "tmp/unique_ior.txt"},
//    {"num_orientations_per_pi", str(num_orientations_per_pi)},
    {"num_z", "2"},
    {"smoothing_distance" , str(smoothing_distance)},
    {"input_table_file", "tmp/contact_table.txt"},
    {"output_table_file", "tmp/table"+str(proc)+".txt"},
  }));
  table->run(mc.get());
  int ior = 0;
  for (int iorall = 0; iorall < 3; ++iorall) {
    if (table->rotator().ior_in_proc(iorall)) {
      double expect1 = -1., expect2 = -1.;
      if (iorall == 0) {
        expect1 = -6.8349628448486328;
        //expect1 = -6.835550;
        expect2 = -0.001538283;
      } else if (iorall == 1) {
        expect1 = -6.7353038787841797;
        //expect1 = -6.736353;
        expect2 = 0.001778595;
      } else if (iorall == 2) {
        expect1 = 0.;
        expect2 = 0.;
      }
      EXPECT_NEAR(expect1, table->rotator().energy_[ior][0], 1e-4);
      EXPECT_NEAR(expect2, table->rotator().energy_[ior][1], 1e-4);
      ++ior;
    }
  }
  if (proc == num_proc - 1) {
    table->combine_table({{"prefix", "tmp/table"},
                          {"suffix", ".txt"}});
  }
}

TEST(MonteCarlo, tabulate_rigid_bodies_LONG) {
  const int num_orientations_per_pi = 1;
  const int num_proc = 2;
  auto mc = MakeMonteCarlo({{
    {"Configuration", {
      {"cubic_side_length", "2e2"},
      {"particle_type0", "../particle/spce.fstprt"},
      {"particle_type1", "../particle/spce.fstprt"},
      //{"particle_type0", "../plugin/aniso/test/data/fc.fstprt"},
      //{"particle_type1", "../plugin/aniso/test/data/fc.fstprt"},
      {"add_particles_of_type0", "1"},
      {"add_particles_of_type1", "1"},
      {"group0", "fixed"},
      {"fixed_particle_type", "0"},
      {"group1", "mobile"},
      {"mobile_particle_type", "1"},
    }},
    {"Potential", {{"Model", "HardSphere"}}}
  }});
  auto table_ior = std::make_shared<TabulateTwoRigidBody3D>(argtype({
    {"num_orientations_per_pi", str(num_orientations_per_pi)},
    {"output_orientation_file", "tmp/unique_ior.txt"}}));
  table_ior->run(mc.get());
  //run(num_orientations_per_pi, 0, 1);
  for (int proc = 0; proc < num_proc; ++proc) {
    run_hs(num_orientations_per_pi, proc, num_proc);
  }
  for (int proc = 0; proc < num_proc; ++proc) {
    run(num_orientations_per_pi, proc, num_proc);
  }
  EXPECT_EQ(table_ior->rotator().fraction_unique(), 1./9.);
  EXPECT_EQ(table_ior->rotator().unique_[0], -1);
  EXPECT_EQ(table_ior->rotator().unique_[1], -1);
  EXPECT_EQ(table_ior->rotator().unique_[2], 0);
  EXPECT_EQ(table_ior->rotator().unique_[3], -1);
  EXPECT_EQ(table_ior->rotator().unique_[4], -1);
  EXPECT_EQ(table_ior->rotator().unique_[5], 3);
}

//TEST(MonteCarlo, analyze_orientations) {
//  auto mc = MakeMonteCarlo({{
//    {"Configuration", {
//      {"cubic_side_length", "2e2"},
//      {"particle_type0", "../particle/spce.fstprt"},
//      {"particle_type1", "../particle/spce.fstprt"},
//      //{"particle_type1", "../particle/propane.fstprt"},
//      //{"particle_type0", "../plugin/aniso/test/data/fc.fstprt"},
//      //{"particle_type1", "../plugin/aniso/test/data/fc.fstprt"},
//      {"add_particles_of_type0", "1"},
//      {"add_particles_of_type1", "1"},
//      {"group0", "fixed"},
//      {"fixed_particle_type", "0"},
//      {"group1", "mobile"},
//      {"mobile_particle_type", "1"},
//    }},
//    {"Potential", {{"Model", "HardSphere"}}}
//  }});
//  auto table_ior = MakeTabulateTwoRigidBody3D({
//    {"input_orientation_file", "../../open_mab_cg/orientations/orientations16.txt"}});
//    //{"input_orientation_file", "../../open_mab_cg/orientations/orientations12_ij.txt"}});
//    //{"input_orientation_file", "tmp/unique_ior.txt"}});
//  table_ior->read_input_orientations_(mc->configuration());
//  INFO(table_ior->rotator().fraction_unique());
//  INFO(table_ior->rotator().num_orientations());
//}

}  // namespace feasst
