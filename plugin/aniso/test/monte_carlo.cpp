#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

TEST(MonteCarlo, VisitModelInnerTable) {
  //const std::string table_file = "../plugin/aniso/tutorial/dat.txt";
  //const std::string table_file = "../plugin/aniso/test/data/dat_3rel_2z.txt";
  const std::string table_file = "../plugin/aniso/test/data/dat_sqw_3rel_2z.txt";
  //const std::string table_file = "../plugin/aniso/test/data/dat_sqw_6rel_2z.txt";
  auto vis = MakeVisitModelInnerTable({{"table_file", table_file}});
//  EXPECT_NEAR(vis->outer()[0][0].minimum(), 1.5, NEAR_ZERO);
//  EXPECT_NEAR(vis->outer()[0][0].maximum(), 1.5, NEAR_ZERO);
  EXPECT_NEAR(vis->inner()[0][0].minimum(), 1.0, NEAR_ZERO);
  EXPECT_NEAR(vis->inner()[0][0].maximum(), 1.0, NEAR_ZERO);
  EXPECT_NEAR(vis->energy()[0][0].minimum(), -1, 1e-3);
  EXPECT_NEAR(vis->energy()[0][0].maximum(), -1., 1e-3);
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../plugin/aniso/particle/aniso_tabular.fstprt"},
      {"xyz_file", "../plugin/aniso/test/data/two.xyz"}}},
    //{"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../plugin/aniso/particle/aniso_tabular.fstprt"},
    //  {"add_particles_of_type0", "1"}}},
    //{"Potential", {{"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    {"Potential", {{"Model", "TwoBodyTable"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    //{"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    //{"Potential", {{"VisitModelInner", "VisitModelInnerTable"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialRotate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"file_name", "tmp/aniso.fst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"file_name", "tmp/aniso.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"file_name", "tmp/aniso.xyze"}, {"euler", "true"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }});
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
                       {"particle_type0", "../plugin/aniso/test/data/rigid_and_connector.fstprt"}}},
    {"Potential", {{"Model", "TwoBodyTable"}, {"VisitModelInner", "VisitModelInnerTable"}, {"table_file", table_file}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1"}, {"tunable_param", "4"}}},
    {"TrialRotate", {{"weight", "1"}, {"tunable_param", "180"}}},
    //{"Run", {{"num_trials", "1"}}},
    {"Run", {{"num_trials", "4"}}},
    {"RemoveTrial", {{"name", "TrialTranslate"}}},
    {"RemoveTrial", {{"name", "TrialRotate"}}},
    {"TrialGrowFile", {{"file_name", "../plugin/aniso/test/data/rigid_and_connector.txt"}}},
    {"Log", {{"trials_per_write", "1"}, {"file_name", "tmp/connector.txt"}}},
    {"Movie", {{"trials_per_write", "1"}, {"file_name", "tmp/connector.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e2)}, {"tolerance", "1e-8"}}},
  }});
  mc->attempt(1e3);
  MonteCarlo mc2 = test_serialize(*mc);
}

}  // namespace feasst
