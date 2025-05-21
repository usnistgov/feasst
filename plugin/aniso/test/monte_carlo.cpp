#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

TEST(MonteCarlo, VisitModelInnerTable) {
  //const std::string table_file = "../plugin/aniso/tutorial/dat.txt";
  //const std::string table_file = "../plugin/aniso/test/data/dat_3rel_2z.txt";
  const std::string table_file = "../plugin/aniso/test/data/dat_sqw_3rel_2z.txt";
  //const std::string table_file = "../plugin/aniso/test/data/dat_sqw_6rel_2z.txt";
  auto vis = std::make_shared<VisitModelInnerTable>(argtype({{"table_file", table_file}}));
//  EXPECT_NEAR(vis->outer()[0][0].minimum(), 1.5, NEAR_ZERO);
//  EXPECT_NEAR(vis->outer()[0][0].maximum(), 1.5, NEAR_ZERO);
  auto config = MakeConfiguration({{"particle_type0", "../particle/atom.txt"}});
  vis->precompute(config.get());
  EXPECT_NEAR(config->table5d()[0][0]->minimum(), 1.0, NEAR_ZERO);
  EXPECT_NEAR(config->table5d()[0][0]->maximum(), 1.0, NEAR_ZERO);
  EXPECT_NEAR(config->table6d()[0][0]->minimum(), -1, 1e-3);
  EXPECT_NEAR(config->table6d()[0][0]->maximum(), -1., 1e-3);
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../plugin/aniso/particle/aniso_tabular.txt"},
      {"xyz_file", "../plugin/aniso/test/data/two.xyz"}}},
    //{"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../plugin/aniso/particle/aniso_tabular.txt"},
    //  {"add_particles_of_type0", "1"}}},
    //{"Potential", {{"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    {"Potential", {{"Model", "TwoBodyTable"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    //{"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    //{"Potential", {{"VisitModelInner", "VisitModelInnerTable"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialRotate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/aniso.fst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/aniso.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/aniso.xyze"}, {"euler", "true"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }}, true);
  EXPECT_EQ(2, mc->configuration().num_particles());
  // EXPECT_NEAR(-1, mc->criteria().current_energy(), NEAR_ZERO);
  EXPECT_TRUE(mc->configuration().unique_type(0).site(0).is_anisotropic());
  EXPECT_FALSE(mc->configuration().particle_type(0).site(0).is_anisotropic());
  EXPECT_TRUE(mc->configuration().particle(0).site(0).is_anisotropic());
//  EXPECT_NEAR(-2.060346185437E+00, mc->configuration().particle(0).site(0).position().coord(0), NEAR_ZERO);
//  EXPECT_TRUE(std::abs(mc->configuration().particle(0).site(0).position().coord(0)-1.077169909511E+00)>1e-8);
}

TEST(MonteCarlo, rigid_body_connector) {
  const std::string table_file = "../plugin/aniso/test/data/dat_sqw_3rel_2z.txt";
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "1672847223"}}},
    //{"RandomMT19937", {{"seed", "time"}}},
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "8"}, {"add_particles_of_type0", "1"},
                       {"particle_type0", "../plugin/aniso/test/data/rigid_and_connector.txt"}}},
    {"Potential", {{"Model", "TwoBodyTable"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1"}, {"tunable_param", "4"}}},
    {"TrialRotate", {{"weight", "1"}, {"tunable_param", "180"}}},
    //{"Run", {{"num_trials", "1"}}},
    {"Run", {{"num_trials", "4"}}},
    {"Remove", {{"name0", "TrialTranslate"}, {"name1", "TrialRotate"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/aniso/test/data/rigid_and_connector_grow.txt"}}},
    {"Log", {{"trials_per_write", "1"}, {"output_file", "tmp/connector.txt"}}},
    {"Movie", {{"trials_per_write", "1"}, {"output_file", "tmp/connector.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e2)}, {"tolerance", "1e-8"}}},
  }}, true);
  mc->attempt(1e3);
  auto mc2 = test_serialize_unique(*mc);
  mc2->attempt(1e3);
}

TEST(MonteCarlo, 4lyt_smoothing) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {
      {"particle_type0", "../plugin/aniso/particle/aniso_tabular.txt"},
      {"xyz_euler_file", "../plugin/aniso/test/data/two.xyze"}}},
    {"Potential", {{"Model", "TwoBodyTable"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", "../plugin/aniso/test/data/table.txt"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
//    {"Log", {{"trials_per_write", "1"}, {"output_file", "tmp/4lyt.txt"}}},
//    {"Run", {{"num_trials", "4"}}},
  }}, true);
//  mc->attempt(1);
  // HWH add serialization test here as well
//  MonteCarlo mc2 = test_serialize(*mc);
//  //mc2.initialize_criteria();
//  EXPECT_NEAR(mc->get_system()->energy(),
//              mc2.get_system()->energy(), NEAR_ZERO);
  EXPECT_NEAR(0.0003444443315056451, mc->get_system()->energy(), NEAR_ZERO);

//  const double temperature = 298.15;
//  const double temp_cel = temperature - 273.15;
//  const double dielectric_water = 87.74 - 0.40008*temp_cel + 9.398e-4*std::pow(temp_cel, 2) - 1.4e-6*std::pow(temp_cel, 3);
//  const double eps_0 = CODATA2018().permitivity_vacuum();
//  const double elem_q = CODATA2018().elementary_charge();
//  const double na = CODATA2018().avogadro_constant();
//  const double kb = CODATA2018().boltzmann_constant();
//  const double ionic_strength = 0.15;
//  const double kappa = std::sqrt(2*std::pow(elem_q,2)*ionic_strength*(1e3)*na/(dielectric_water*eps_0*kb*temperature*1e20));
//  const double cutoff = 5./kappa;
//  auto all_atom = MakeMonteCarlo({{
//    {"Configuration", {
//      {"cubic_side_length", "500"},
//      {"cutoff", str(cutoff)},
//      {"particle_type0", "../../open_mab_cg/4lyt_ph6/4lyt.txt"},
//      {"add_particles_of_type0", "2"}}},
//    {"Potential", {{"Model", "ModelTwoBodyFactory"}, {"model0", "HardSphere"}, {"model1", "LennardJones"}, {"model2", "DebyeHuckel"}, {"kappa", str(kappa)}, {"dielectric", str(dielectric_water)}, {"smoothing_distance", "2"}, {"energy_cutoff", "1e100"}}},
//    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1."}}},
//    {"Metropolis", {{}}},
////    {"Log", {{"trials_per_write", "1"}, {"output_file", "tmp/4lyt.txt"}}},
////    {"Run", {{"num_trials", "4"}}},
//  }}, true);
//  Select sel;
//  sel.add_particle(mc->configuration().particle(1), 1);
//  Position disp;
//  //disp.set_vector({30, 0, 0});
//  //disp.set_vector({2.870980e+01, 0, 0});
//  disp.set_vector({65.933407323399, 0, 0});
//  //disp.set_vector({67.933407323399, 0, 0});
//  all_atom->get_system()->get_configuration()->displace_particles(sel, disp);
//  all_atom->initialize_criteria();
//  EXPECT_NEAR(0, all_atom->criteria().current_energy(), NEAR_ZERO);
}

}  // namespace feasst
