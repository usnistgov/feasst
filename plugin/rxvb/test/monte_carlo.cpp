#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "actions/include/remove.h"
#include "actions/include/run.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/wang_landau.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "morph/include/trial_morph.h"
#include "morph/include/trial_morph_expanded.h"
#include "morph/include/macrostate_morph.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"

namespace feasst {

// pure ideal gas starting with 0 and 1 should give <N_2>=0.5

// Testing:
// 2 particle (0&1) ideal gas: Morph and RXNAVB both give <N> = 0.5
// Add SquareWell and Confinement with -1 mu3: Morph and RXNAVB both give ~0.33
// Add a third particle? Morph ~ 0.47, RXNAVB ~ 0.375
// Remove confinement slab: Morph ~ 0.47, RXNAVB ~ 0.375
// Remove translates: Morph ~ 0.425 , RXNAVB ~ 0.425
// Use 3-sep without translates: Morph~ 0.53, RXNAVB ~ 0.5 .. NOT a fair test because Morph can rx with distant particle
// Added both Morph and RXNAVB: ~0.5256 +/- 0.0002, 0.52586 +/- 0.00020, Morph alone: 0.5317 +/- 0.0002
// Added translates back, then comparing both Morph and RXNAVB: 0.4627 pm 0.001, Morph only: 0.46649 pm 0.0009
// 1e7->1e8 trials: Both: 0.465424290 pm 0.00036 (z~1.7), Morph only: 0.46647213 pm 0.00047, RXN only: 0.37865451999999999 pm 0.002487
// Reduce Morph weights to 0.1, Both: 0.459658459 pm 0.00023 (z~13)
// Test to see if Morph/RXNAVB differences disappear if there are no fixed sites.
// Added translates for 0/2 and dropped trials to 1e7: Morph: 0.468585699 pm 0.00109, RXNAVB: 0.386549 pm 0.00618
// translates weight_per_num_frac ->weight,RXNAVB: 0.38534, Morph 0.46876
// Rederived symmetric
// RXNAVB: 0.47174 pm 0.01, Morph: 0.46880 pm 0.00132

TEST(MonteCarlo, rxvb_ideal) {
  const std::string tpis("1e2");
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "1234"}}},
    {"Configuration", {{"particle_type", "../particle/atom.txt,../particle/atom.txt,../particle/atom.txt,../particle/atom.txt"},
                       {"sigma", "0"},
                       {"epsilon", "0"}, {"epsilon0_1", "1"}, {"epsilon2_3", "2"},
                       {"cutoff", "0"}, {"cutoff0_1", "1.5"}, {"cutoff2_3", "1.5"},
                       {"xyz_file", "../plugin/rxvb/test/data/chem2.xyz"},
                       //{"xyz_file", "../plugin/rxvb/test/data/chem3_sep.xyz"},
                       {"group0", "fluid"}, {"fluid_particle_type", "1"}, {"fluid_particle_type", "3"},
                      }},
    {"NeighborCriteria", {{"site_type0", "0"}, {"site_type1", "1"}, {"site_type0_alt", "2"}, {"site_type1_alt", "3"},
                          {"maximum_distance", "1.5"}, {"minimum_distance", "1"}}},
    //{"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapAllCriteria"}}},
    //{"Potential", {{"Model", "SquareWell"}}},
    //{"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}}},
    {"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}, {"VisitModel", "VisitModelCell"}, {"min_length", "max_cutoff"}}},
    //{"Potential", {{"Model", "IdealGas"}, {"EnergyMap", "EnergyMapAllCriteria"}}},
    //{"Potential", {{"Model", "IdealGas"}}},
    //{"Potential", {{"Model", "IdealGas"}, {"EnergyMap", "EnergyMapNeighbor"}}},
    //{"Potential", {{"Model", "IdealGas"}, {"EnergyMap", "EnergyMapNeighborCriteria"}}},
    //{"Potential", {{"Model", "ModelHardShape"}, {"shape_file", "../plugin/rxvb/test/data/chem_shape_file.txt"}, {"group", "fluid"}}},
    {"ThermoParams", {{"beta", "1e-15"}, {"chemical_potential", "0.,0.,0.,0."}}},
    {"Metropolis", {{}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "3"}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "0"}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "2"}}},
    {"TrialTranslate", {{"weight", "1."}, {"particle_type", "1"}}},
    {"TrialTranslate", {{"weight", "1."}, {"particle_type", "3"}}},
    {"TrialTranslate", {{"weight", "1."}, {"particle_type", "0"}}},
    {"TrialTranslate", {{"weight", "1."}, {"particle_type", "2"}}},
    {"TrialMorph", {{"weight", "0.2"}, {"particle_type", "0,1"}, {"particle_type_morph", "2,3"}}},
//    {"TrialRxVB", {{"target_particle_type", "0"}, {"target_particle_type_morph", "2"}, {"particle_type", "1"}, {"particle_type_morph", "3"}}},
    {"CheckEnergy", {{"trials_per_update", tpis}, {"decimal_places", "8"}}},
    {"NumParticles", {{"trials_per_write", tpis}, {"output_file", "tmp/chemn.csv"}, {"particle_type", "2"}}},
    {"Energy", {{"trials_per_write", tpis}, {"output_file", "tmp/chem_en.csv"}}},
    {"Log", {{"trials_per_write", tpis}, {"output_file", "tmp/chem.csv"}}},
    {"Movie", {{"trials_per_write", tpis}, {"output_file", "tmp/chem.xyz"}}},
    //{"ProfileCPU", {{"trials_per_write", tpis}, {"output_file", "tmp/ljp.csv"}, {"trials_per_update", "1e3"}}},
  }}, false); //true);
  for (int trial = 0; trial < 1e4; ++trial) {
    mc->run_num_trials(1);
//    const double en = mc->criteria().current_energy();
//    if (en > -1) {
//      DEBUG("en " << en);
//      const Particle& part0 = mc->configuration().particle(0);
//      DEBUG("part0 xyz " << part0.site(0).position().str());
//      DEBUG("part0 type " << part0.type());
//      const Particle& part1 = mc->configuration().particle(1);
//      DEBUG("part1 xyz " << part1.site(0).position().str());
//      DEBUG("part1 type " << part1.type());
//    }
  }
  const Accumulator& acc = mc->analyze(0).accumulator();
  INFO(acc.str());
  EXPECT_NEAR(acc.average(), 0.5, 5*acc.block_stdev());
}

std::shared_ptr<MonteCarlo> make_mc(const std::string xyz, const double mu, const std::string fstprt = "../particle/atom.txt") {
  const std::string tpis("1");
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "1234"}}},
    {"Configuration", {{"particle_type", fstprt+','+fstprt+','+fstprt+','+fstprt},
                       {"epsilon", "0"}, {"epsilon0_1", "1"}, {"epsilon2_3", "2"},
                       {"cutoff", "1"}, {"cutoff0_1", "1.5"}, {"cutoff2_3", "1.5"},
                       {"xyz_file", xyz},
                       {"group0", "fluid"}, {"fluid_particle_type0", "1"}, {"fluid_particle_type1", "3"},
                       }},
    {"WriteModelParams", {{"output_file", "chem_model_params.txt"}}},
    {"NeighborCriteria", {{"site_type0", "0"}, {"site_type1", "1"}, {"site_type0_alt", "2"}, {"site_type1_alt", "3"},
                          {"maximum_distance", "1.5"}, {"minimum_distance", "1"}}},
    {"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}, {"VisitModel", "VisitModelCell"}, {"min_length", "max_cutoff"}}},
    //{"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}}},
    {"Potential", {{"Model", "ModelHardShape"}, {"shape_file", "../plugin/rxvb/test/data/chem_shape_file.txt"}, {"group", "fluid"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "0.,0.,0.,"+str(mu)}}},
    {"Metropolis", {{}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "3"}}},
//    {"TrialAdd", {{"particle_type", "1"}}},
//    {"Run", {{"until_num_particles", "10"}, {"particle_type", "1"}}},
//    {"Remove", {{"name", "TrialAdd"}}},
//    {"TrialMorph", {{"weight", "1"},
//                    {"particle_type0", "0"}, {"particle_type_morph0", "2"},
//                    {"particle_type1", "1"}, {"particle_type_morph1", "3"}}},
//    {"TrialMorph", {{"weight", "1"},
//                    {"particle_type0", "2"}, {"particle_type_morph0", "0"},
//                    {"particle_type1", "3"}, {"particle_type_morph1", "1"}}},
    {"TrialRxVB", {{"target_particle_type", "0"}, {"target_particle_type_morph", "2"}, {"particle_type", "1"}, {"particle_type_morph", "3"}}},
    {"CheckEnergy", {{"trials_per_update", tpis}, {"decimal_places", "8"}}},
    //{"Log", {{"trials_per_write", tpis}, {"output_file", "tmp/lj.csv"}}},
    //{"Movie", {{"trials_per_write", tpis}, {"output_file", "tmp/lj.xyz"}}},
    {"NumParticles", {{"trials_per_write", tpis}, {"output_file", "tmp/chemn.csv"}, {"particle_type", "2"}}},
    //{"ProfileCPU", {{"trials_per_write", tpis}, {"output_file", "tmp/ljp.csv"}, {"trials_per_update", "1e3"}}},
  }}, true);
  return mc;
}

TEST(MonteCarlo, rxvb2) {
  std::shared_ptr<MonteCarlo> mc = make_mc("../plugin/rxvb/test/data/chem2.xyz", -1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 0);
  mc->run_num_trials(1);
  // by seed, the first reaction is AVB and accepted
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 1);
  mc->run_num_trials(1);
  // by seed, the second reaction is noAVB and accepted
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 0);
}

TEST(MonteCarlo, rxvb3) {
  std::shared_ptr<MonteCarlo> mc = make_mc("../plugin/rxvb/test/data/chem3.xyz", -0.5);
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 2);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 0);
  mc->run_num_trials(1);
  // by seed, the first reaction is AVB and accepted
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 1);
  mc->run_num_trials(1);
  // by seed, the second reaction is noAVB and accepted
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 2);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 0);
}

TEST(MonteCarlo, rxvb3_sep) {
  std::shared_ptr<MonteCarlo> mc = make_mc("../plugin/rxvb/test/data/chem3_sep.xyz", -0.5);
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 2);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 0);
  mc->run_num_trials(1);
  // by seed, the first reaction is AVB and accepted
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 1);
  mc->run_num_trials(1);
  // by seed, the second reaction is noAVB and accepted
  EXPECT_EQ(mc->configuration().num_particles_of_type(0), 0);
  EXPECT_EQ(mc->configuration().num_particles_of_type(1), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(2), 1);
  EXPECT_EQ(mc->configuration().num_particles_of_type(3), 1);
}

std::shared_ptr<MonteCarlo> make_mc_dimer(const std::string xyz, const double mu, const std::string fstprt = "../plugin/rxvb/particle/dimer.txt") {
  const std::string tpis("1");
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "1234"}}},
    {"Configuration", {{"particle_type0", fstprt},
                       {"particle_type1", fstprt},
                       {"particle_type2", fstprt},
                       {"particle_type3", fstprt},
                       {"epsilon", "1"}, {"epsilon4_6", "2"},
                       {"cutoff", "1.5"},
                       {"xyz_file", xyz},
                       {"group0", "fluid"}, {"fluid_particle_type0", "1"}, {"fluid_particle_type1", "3"},
                       }},
    {"WriteModelParams", {{"output_file", "chem_model_params.txt"}}},
    {"NeighborCriteria", {{"site_type0", "0"}, {"site_type1", "2"}, {"site_type0_alt", "4"}, {"site_type1_alt", "6"},
                          {"maximum_distance", "1.5"}, {"minimum_distance", "1"}}},
    {"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}, {"VisitModel", "VisitModelCell"}, {"min_length", "max_cutoff"}}},
    //{"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}}},
    {"Potential", {{"Model", "ModelHardShape"}, {"shape_file", "../plugin/rxvb/test/data/chem_shape_file.txt"}, {"group", "fluid"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "0."}, {"chemical_potential1", "0."},
                                     {"chemical_potential2", "0."}, {"chemical_potential3", str(mu)}}},
    {"Metropolis", {{}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}}},
//    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "3"}}},
//    {"TrialAdd", {{"particle_type", "1"}}},
//    {"Run", {{"until_num_particles", "10"}, {"particle_type", "1"}}},
//    {"Remove", {{"name", "TrialAdd"}}},
//    {"TrialMorph", {{"weight", "1"},
//                    {"particle_type0", "0"}, {"particle_type_morph0", "2"},
//                    {"particle_type1", "1"}, {"particle_type_morph1", "3"}}},
//    {"TrialMorph", {{"weight", "1"},
//                    {"particle_type0", "2"}, {"particle_type_morph0", "0"},
//                    {"particle_type1", "3"}, {"particle_type_morph1", "1"}}},
    {"TrialRxVB", {{"target_particle_type", "0"}, {"target_particle_type_morph", "2"}, {"particle_type", "1"}, {"particle_type_morph", "3"}, {"reference_index", "0"}}},
    {"CheckEnergy", {{"trials_per_update", tpis}, {"decimal_places", "8"}}},
    //{"Log", {{"trials_per_write", tpis}, {"output_file", "tmp/lj.csv"}}},
    //{"Movie", {{"trials_per_write", tpis}, {"output_file", "tmp/lj.xyz"}}},
    {"NumParticles", {{"trials_per_write", tpis}, {"output_file", "tmp/chemn.csv"}, {"particle_type", "2"}}},
    //{"ProfileCPU", {{"trials_per_write", tpis}, {"output_file", "tmp/ljp.csv"}, {"trials_per_update", "1e3"}}},
  }}, true);
  return mc;
}
TEST(MonteCarlo, rxvb2_dimer) {
  std::shared_ptr<MonteCarlo> mc = make_mc_dimer("../plugin/rxvb/test/data/chem2_dimer.xyz", -1);
  mc->get_random()->seed(1234567);
  const Configuration& config = mc->configuration();
  Position p0s0 = config.particle(0).site(0).position();
  Position p0s1 = config.particle(0).site(1).position();
  Position p1s0 = config.particle(1).site(0).position();
  Position p1s1 = config.particle(1).site(1).position();
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(0).site(1).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(1).site(1).position().str());
  EXPECT_EQ(config.num_particles_of_type(0), 1);
  EXPECT_EQ(config.num_particles_of_type(1), 1);
  EXPECT_EQ(config.num_particles_of_type(2), 0);
  EXPECT_EQ(config.num_particles_of_type(3), 0);
  mc->run_num_trials(1);
  // by seed, the first reaction is AVB and accepted
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(0).site(1).position().str());
  DEBUG(config.particle(1).site(0).position().str());
  DEBUG(config.particle(1).site(1).position().str());
  EXPECT_TRUE(config.particle(0).site(0).position().is_equal(p0s0));
  EXPECT_FALSE(config.particle(0).site(1).position().is_equal(p0s1));
  EXPECT_TRUE(config.particle(1).site(0).position().is_equal(p1s0));
  EXPECT_FALSE(config.particle(1).site(1).position().is_equal(p1s1));
  EXPECT_EQ(config.num_particles_of_type(0), 0);
  EXPECT_EQ(config.num_particles_of_type(1), 0);
  EXPECT_EQ(config.num_particles_of_type(2), 1);
  EXPECT_EQ(config.num_particles_of_type(3), 1);
}

TEST(MonteCarlo, rxvb_dimer) {
  //const std::string tpis("1e0");
  const std::string tpis("1e2");
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"particle_type0", "../plugin/rxvb/particle/dimer.txt"},
                       {"particle_type1", "../plugin/rxvb/particle/dimer.txt"},
                       {"particle_type2", "../plugin/rxvb/particle/dimer.txt"},
                       {"particle_type3", "../plugin/rxvb/particle/dimer.txt"},
                       {"epsilon", "1"}, {"epsilon4_6", "2"},
                       {"cutoff", "1.5"},
                       {"xyz_file", "../plugin/rxvb/test/data/chem2_dimer.xyz"},
                       {"group0", "fluid"}, {"fluid_particle_type0", "1"}, {"fluid_particle_type1", "3"},
                      }},
    {"NeighborCriteria", {{"site_type0", "0"}, {"site_type1", "2"}, {"site_type0_alt", "4"}, {"site_type1_alt", "6"},
                          {"maximum_distance", "1.5"}, {"minimum_distance", "1"}}},
    {"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}}},
    //{"Potential", {{"Model", "SquareWell"}, {"EnergyMap", "EnergyMapNeighborCriteria"}, {"VisitModel", "VisitModelCell"}, {"min_length", "max_cutoff"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "0."}, {"chemical_potential1", "0."},
                                     {"chemical_potential2", "0."}, {"chemical_potential3", "0."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "3"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "0"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "2"}}},
//    {"TrialMorph", {{"weight", "0.1"},
//                    {"particle_type0", "0"}, {"particle_type_morph0", "2"},
//                    {"particle_type1", "1"}, {"particle_type_morph1", "3"}}},
//    {"TrialMorph", {{"weight", "0.1"},
//                    {"particle_type0", "2"}, {"particle_type_morph0", "0"},
//                    {"particle_type1", "3"}, {"particle_type_morph1", "1"}}},
    {"TrialAdd", {{"particle_type", "1"}}},
    {"Run", {{"until_num_particles", "20"}, {"particle_type", "1"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"TrialRxVB", {{"target_particle_type", "0"}, {"target_particle_type_morph", "2"}, {"particle_type", "1"}, {"particle_type_morph", "3"}, {"reference_index", "0"}}},
    {"CheckEnergy", {{"trials_per_update", tpis}, {"decimal_places", "8"}}},
    {"NumParticles", {{"trials_per_write", tpis}, {"output_file", "tmp/chemn.csv"}, {"particle_type", "2"}}},
    {"Energy", {{"trials_per_write", tpis}, {"output_file", "tmp/chem_en.csv"}}},
    {"Log", {{"trials_per_write", tpis}, {"output_file", "tmp/chem.csv"}}},
    {"Movie", {{"trials_per_write", tpis}, {"output_file", "tmp/chem.xyz"}}},
    //{"ProfileCPU", {{"trials_per_write", tpis}, {"output_file", "tmp/ljp.csv"}, {"trials_per_update", "1e3"}}},
  }}, true);
  for (int trial = 0; trial < 1e2; ++trial) {
    mc->run_num_trials(1);
  }
  //const Accumulator& acc = mc->analyze(0).accumulator();
  //INFO(acc.str());
  //EXPECT_NEAR(acc.average(), 0.5, 5*acc.block_stdev());
}

TEST(MonteCarlo, TrialSwapPosition) {
  const std::string tpc = "1e1";
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "20"}, {"add_particles_of_type0", "1"},
        {"particle_type0", "../particle/dimer.txt"}, {"particle_type1", "../particle/spce.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "10"}, {"chemical_potential1", "10"}}},
    {"Metropolis", {{}}},
    {"TrialAdd", {{"particle_type", "1"}}},
    {"Run", {{"until_num_particles", "3"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Metropolis", {{"trials_per_cycle", tpc}, {"cycles_to_complete", "1e2"}}},
    {"TrialSwapPosition", {{"particle_types", "0,1"}}},
    {"Tune", {{"trials_per_tune", "4"}}},
    {"Log", {{"trials_per_write", tpc}, {"output_file", "tmp/swap.csv"}}},
    {"Movie", {{"trials_per_write", tpc}, {"output_file", "tmp/swapc.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", tpc}, {"tolerance", str(1e-9)}}},
    {"Run", {{"until", "complete"}}},
  }}, true);
  auto mc2 = test_serialize_unique(*mc);
  mc2->run_num_trials(1e2);
}

}  // namespace feasst
