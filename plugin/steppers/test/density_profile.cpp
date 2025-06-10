#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/density_profile.h"

namespace feasst {

TEST(DensityProfile, serialize) {
  auto density_profile = MakeDensityProfile({{"output_file", "tmp"}});
  auto density_profile2 = test_serialize<DensityProfile, Analyze>(*density_profile);
}

TEST(DensityProfile, dimer) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "4.5"}, {"particle_type0", "../particle/dimer.txt"}, {"cutoff", "1"}}},
    //{"Configuration", {{"cubic_side_length", "4"}, {"particle_type0", "../particle/heterodimer.txt"}, {"cutoff", "1"}}},
    {"Potential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "3"}}},
    {"TrialParticlePivot", {{"particle_type", "0"}, {"tunable_param", "3"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "10"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Log", {{"trials_per_write", "100"}, {"output_file", "tmp/dimer.csv"}}},
    {"Movie", {{"trials_per_write", "100"}, {"output_file", "tmp/dimer.xyz"}}},
    }});
  auto profile = MakeDensityProfile({{"trials_per_update", "100"}, {"dr", "0.5"},
    {"trials_per_write", "1000"}, {"dimension", "2"}, {"output_file", "tmp/dimerprof.txt"}});
  mc->add(profile);
  mc->attempt(1e4);
  auto profile2 = test_serialize(*profile);
  const int num = static_cast<int>(profile2.profile().size());
  EXPECT_EQ(num, 11);
  EXPECT_NEAR(profile2.profile()[0][0][0], -2.5, NEAR_ZERO);
  EXPECT_NEAR(profile2.profile()[0][0][1], 0., NEAR_ZERO);
  EXPECT_NEAR(profile2.profile()[1][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[2][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[3][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[4][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[5][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[6][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[7][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[8][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[9][0][1], 1./(num-2), 0.035);
  EXPECT_NEAR(profile2.profile()[10][0][1], 0., NEAR_ZERO);
  EXPECT_NEAR(profile2.profile()[10][0][0], 2.5, NEAR_ZERO);
}

}  // namespace feasst
