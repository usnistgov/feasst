#include "utils/test/utils.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "actions/include/record.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

TEST(Record, serialize) {
  auto obj = std::make_shared<Record>(argtype({{"save_positions", "something"}}));
  Action obj2 = test_serialize(*obj);
}

TEST(Record, save_load) {
  std::string seed = "-1";
  //std::string seed = "123";
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", seed}}},
    {"Configuration", {{"side_length", "15,15"},
      {"particle_type", "pt:../plugin/aniso/particle/aniso_tabular2d3l.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "1"}}},
    {"TrialRotate", {{"particle_type", "pt"}, {"tunable_param", "100"}}},
    {"TrialAdd", {{"particle_type", "pt"}}},
    {"Run", {{"until_num_particles", "2"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Record", {{"save_positions", "tmp/rec.xyze"}}},
  }}, true);
  auto mc2 = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", seed}}},
    {"Configuration", {{"xyz_euler_file", "tmp/rec.xyze"},
      {"particle_type", "pt:../plugin/aniso/particle/aniso_tabular2d3l.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Record", {{"save_positions", "tmp/rec2.xyze"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1"}}},
    {"Metropolis", {{}}},
  }}, true);
  for (int part = 0; part < 2; ++part) {
    for (int site = 0; site < 2; ++site) {
      const Site& s1 =  mc->configuration().particle(part).site(site);
      const Site& s2 = mc2->configuration().particle(part).site(site);
      for (int dim = 0; dim < mc->configuration().dimension(); ++dim) {
        EXPECT_NEAR(s1.position(dim), s2.position(dim), 1e-4);
      }
      EXPECT_NEAR(s1.euler().phi(), s2.euler().phi(), 1e-4);
      EXPECT_NEAR(s1.euler().theta(), s2.euler().theta(), 1e-4);
      EXPECT_NEAR(s1.euler().psi(), s2.euler().psi(), 1e-4);
    }
  }
}

}  // namespace feasst
