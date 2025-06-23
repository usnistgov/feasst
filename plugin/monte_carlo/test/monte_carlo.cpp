#include <memory>
#include <cmath>
#include "monte_carlo/test/monte_carlo_utils.h"
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "utils/include/checkpoint.h"
#include "utils/include/progress_report.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/dont_visit_model.h"
#include "system/include/ideal_gas.h"
#include "system/include/thermo_params.h"
#include "system/include/potential.h"
#include "system/include/visit_model_inner.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_volume.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/remove.h"
#include "monte_carlo/include/convert_to_ref_potential.h"
#include "monte_carlo/test/monte_carlo_benchmark.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/movie.h"
#include "steppers/include/tune.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/profile_trials.h"
#include "steppers/include/volume.h"
#include "opt_lj/include/visit_model_opt_lj.h"
//#include "cluster/include/energy_map_neighbor.h"
//#include "cluster/include/energy_map_all.h"

namespace feasst {

TEST(MonteCarlo, serialize) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/ljrst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"ProfileCPU", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj_prof.csv"}}},
    {"Tune", {{}}},
  }}, true);
  auto mc2 = test_serialize_unique(*mc);
  mc2->run_num_trials(10);
  EXPECT_EQ(mc2->analyze(0).class_name(), "Log");
  EXPECT_EQ(mc2->analyze(1).class_name(), "Movie");
  EXPECT_EQ(mc2->modify(0).class_name(), "CheckEnergy");
  EXPECT_EQ(mc2->modify(1).class_name(), "Tune");

  auto mc3 = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialTransfer", {{"weight", "4."}, {"particle_type", "0"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/ljrst"}}},
    {"Log", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
  }}, true);

  auto mc4 = test_serialize_unique(*mc3);
  EXPECT_EQ(mc4->trial(0).class_name(), "TrialTranslate");
  EXPECT_EQ(mc4->trial(0).weight(), 1.);
  EXPECT_EQ(mc4->trial(1).class_name(), "TrialAdd");
  EXPECT_EQ(mc4->trial(1).weight(), 2.);
  EXPECT_EQ(mc4->trial(2).class_name(), "TrialRemove");
  EXPECT_EQ(mc4->trial(2).weight(), 2.);
  EXPECT_EQ(mc4->analyze(0).class_name(), "Log");
  EXPECT_EQ(mc4->analyze(1).class_name(), "Movie");
  EXPECT_EQ(mc4->modify(0).class_name(), "CheckEnergy");
  EXPECT_EQ(mc4->modify(1).class_name(), "Tune");
}

TEST(MonteCarlo, NVT_NO_FEASST_BENCHMARK_LONG) {
  srand(123);
  RandomBenchmark random;
  serial(&random);
}

TEST(MonteCarlo, NVT_BENCHMARK_LONG) {
  auto mc = MakeMonteCarlo({{
      {"RandomMT19937", {{"seed", "123"}}},
      {"Configuration", {{"cubic_side_length", "8"},
        {"particle_type0", "../particle/lj.txt"},
        {"xyz_file", "../plugin/monte_carlo/test/data/bench.xyz"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      //{"Potential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelOptLJ"}}},
      {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
  }}, true);
  mc->attempt(1e6);
  //mc.attempt(1e6);  // 2.5s with 50 i9 13900K 7/:1/2024
  //mc.attempt(1e6);  // 4.5s with 50 (see opt_lj for 3s)
  //mc.attempt(1e6);  // 5.4s with 50 (see opt_lj for 4.3s)
  // mc.seek_num_particles(450);
  // mc.attempt(1e5);  // 15 sec with 450 on slow computer
  //mc.attempt(1e3);
  // DEBUG("\n" << mc.timer_str());
}

TEST(MonteCarlo, NVT_cells_BENCHMARK_LONG) {
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "1346867550"}}},
    {"Configuration", {{"cubic_side_length", "12"},
                       {"particle_type", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelCell"}, {"min_length", "3"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {
      {"weight", "1"},
      {"tunable_param", "1"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "200"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
  }}, true);
  mc->attempt(1e5);
}

TEST(MonteCarlo, NVT_cells2_BENCHMARK_LONG) {
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "1346867550"}}},
    {"Configuration", {{"cubic_side_length", "12"},
                       {"particle_type", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {
      {"weight", "1"},
      {"tunable_param", "1"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "200"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
    {"Log", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/cell.txt"}}},
    {"Movie", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/cell.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"OptimizedPotential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelCell"}, {"min_length", "3"}}}
  }}, true);
  mc->initialize_system(0);
  mc->attempt(1e5);
}

TEST(MonteCarlo, NVT_SRSW) {
  const int nMol = 500;
  const double rho = 1e-3;
  const double length = std::pow(static_cast<double>(nMol)/rho, 1./3.);
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", str(length)}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", str(nMol)}}},
    {"Remove", {{"name_contains", "Add"}}},
    {"Log", {{"trials_per_write", "1e3"}, {"output_file", "tmp/lj.csv"}}},
    {"Movie", {{"trials_per_write", "1e3"}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", "1e3"}, {"decimal_places", "9"}}},
    {"Tune", {{}}},
  }}, true);
  Accumulator pe;
  const int num_trials = 1e3;
  for (int trial = 0; trial < num_trials; ++trial) {
    mc->attempt(1);  // ~4 seconds
    pe.accumulate(mc->criteria().current_energy());
  }
  // HWH temperature not set
  DEBUG("pe " << pe.average());
  EXPECT_LE(mc->modify(0).accumulator().average(), 1e-13);
}

TEST(MonteCarlo, GCMC) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "-3"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialTransfer", {{"weight", "4."}, {"particle_type", "0"}}},
    {"NumParticles", {{"trials_per_write", str(1e5)},
                      {"output_file", "tmp/ljnum.txt"}}},
    {"Log", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"ProfileCPU", {{"trials_per_write", str(1e3)},
      {"append", "true"},
      {"output_file", "tmp/lj_profile.txt"}}},
  }}, true);
  const int trials = 1e4;
  //const int trials = 1e6;
  EXPECT_NEAR(mc->trial(0).weight(), 1, NEAR_ZERO);
  EXPECT_NEAR(mc->trial(1).weight(), 2, NEAR_ZERO);
  EXPECT_NEAR(mc->trial(2).weight(), 2, NEAR_ZERO);
  mc->attempt(trials);
  EXPECT_EQ(mc->trials().num(), 3);
  EXPECT_NEAR(mc->trial(0).num_attempts(), trials/5, trials*0.02);
  EXPECT_NEAR(mc->trial(1).num_attempts(), trials*2/5., trials*0.025);
  EXPECT_NEAR(mc->trial(2).num_attempts(), trials*2/5., trials*0.025);
}

TEST(MonteCarlo, GCMC_cell) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"RefPotential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelCell"}, {"min_length", "1"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"CheckEnergy", {{"trials_per_update", "1e4"}, {"tolerance", "1e-9"}}},
    {"Tune", {{}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "-6"}}},
    {"TrialTransfer", {{"particle_type", "0"}, {"num_steps", "4"}, {"reference_index", "0"}}},
    {"NumParticles", {{"trials_per_write", "1e5"},
                      {"output_file", "tmp/ljnum.txt"}}},
  }}, true);
  EXPECT_EQ(mc->system().num_references(), 1);
  //mc->attempt(1e4);
  for (int i = 0; i < 1e4; ++i) {
    mc->attempt(1);
    mc->system().potential(0).visit_model().check(mc->configuration());
  }
}

TEST(MonteCarlo, ConstrainNumParticles) {
  for (const double minimum : {0, 1}) {
    auto mc = MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}, {"add_particles_of_type0", "1"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
      {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
//      {"LogAndMovie", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj"}}},
      {"CheckEnergy", {{"trials_per_update", "1e4"}, {"tolerance", "1e-9"}}},
      {"Tune", {{}}},
      {"ThermoParams", {{"beta", "0.2"}, {"chemical_potential", "-20."}}},
      {"Metropolis", {{"Constraint", "ConstrainNumParticles"},
        {"minimum", str(minimum)}, {"maximum", str(minimum+1)}}},
      {"TrialAddRemove", {{"particle_type", "0"}}},
      {"NumParticles", {{"trials_per_write", "10000"}, {"output_file", "tmp/lj.csv"}}},
    }}, true);
    const int index = mc->num_analyzers() - 1;
    DEBUG("index " << index);
    mc->attempt(14);
    if (minimum == 0) {
      EXPECT_LE(mc->analyze(index).accumulator().average(), 1);
    } else if (minimum == 1) {
      EXPECT_GE(mc->analyze(index).accumulator().average(), 1);
    }
  }
}

TEST(MonteCarlo, ideal_gas_pressure_LONG) {
  const int num = 10;
  const double beta = 1.2, pressure = 0.2;
  const double volume = num/beta/pressure;
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"},
      {"particle_type", "../particle/atom.txt"},
      {"add_particles_of_type0", str(num)},
      {"cutoff", "0"}}},
    {"Potential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", str(beta)}, {"pressure", str(pressure)}}},
    {"Metropolis", {{}}},
    {"TrialVolume", {{"tunable_param", "0.5"}}},
    {"Tune", {{}}},
    {"Volume", {{"trials_per_write", "1e2"},
                {"output_file", "tmp/ideal_gas_volume"}}},
  }}, true);
  mc->attempt(1e6);
  const Accumulator& vol = mc->analyzers()[mc->num_analyzers()-1]->accumulator();
  INFO("volume: " << vol.average() << " +/- " << vol.block_stdev());
  EXPECT_NEAR(vol.average(), volume, 4*vol.block_stdev());
}

TEST(MonteCarlo, lj_npt) {
  const int num = 10;
  const double beta = 1.5, pressure = 0.002;
  const std::string trials_per = "1e0";
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/atom.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", str(beta)}, {"pressure", str(pressure)}, {"chemical_potential", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{}}},
    {"CheckEnergy", {{"trials_per_update", trials_per}, {"tolerance", "1e-4"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", str(num)}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"TrialVolume", {{"tunable_param", "0.5"}}},
    {"Tune", {{}}},
    {"Log", {{"output_file", "tmp/lj.csv"}}},
    {"Volume", {{"trials_per_write", trials_per},
                     {"output_file", "tmp/lj_npt_vol"}}},
  }}, true);
  mc->attempt(1e2);
  const Accumulator& vol = mc->analyzers()[mc->num_analyzers()-1]->accumulator();
  DEBUG("volume: " << vol.average() << " +/- " << vol.block_stdev());
}

TEST(MonteCarlo, arglist_unrecognized) {
  TRY(
    MakeMonteCarlo({{{"Banana", {{}}}}}, true);
    CATCH_PHRASE("Unrecognized argument: Banana");
  );
  TRY(
    MakeMonteCarlo({{{"Metropolis", {{}}}}}, true);
    CATCH_PHRASE("set System before Criteria");
  );
  TRY(
    MakeMonteCarlo({{{"RandomMT19937", {{"this_is_not", "an_expected_argument"}}}}}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{{"Checkpoint", {{"this_is_not", "an_expected_argument"}}}}}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{{"Configuration", {{"particle_type0", "../particle/lj.txt"}, {"this_is_not", "an_expected_argument"}}}}}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"particle_type0", "../particle/lj.txt"}}},
      {"Potential", {{"this_is_not", "an_expected_argument"}}}
    }}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"this_is_not", "an_expected_argument"}}},
    }}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{"this_is_not", "an_expected_argument"}}},
    }}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{"this_is_not", "an_expected_argument"}}},
    }}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{}}},
      {"Energy", {{"this_is_not", "an_expected_argument"}}},
    }}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{}}},
      {"Energy", {{"output_file", "tmp/lj_en.csv"}}},
      {"Tune", {{"this_is_not", "an_expected_argument"}}},
    }}, true);
    CATCH_PHRASE("unused argument");
  );
  TRY(
    MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{}}},
      {"Energy", {{"output_file", "lj_en.csv"}}},
      {"Tune", {{}}},
      {"Run", {{"this_is_not", "an_expected_argument"}}},
    }}, true);
    CATCH_PHRASE("unused argument");
  );
}

TEST(MonteCarlo, argslist_order) {
  auto mc = MakeMonteCarlo({{
    {"RandomModulo", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "8"},
                       {"particle_type0", "../particle/lj.txt"},
                       {"particle_type1", "../particle/atom.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "0.1"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
  }}, true);
  EXPECT_EQ(1.2, mc->thermo_params().beta());
}

TEST(MonteCarlo, arglist) {
  auto mc = MakeMonteCarlo({{
    {"Checkpoint", {{"checkpoint_file", "tmp/lj.fst"}}},
    {"RandomModulo", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "8"},
                       {"particle_type0", "../particle/lj.txt"},
                       {"particle_type1", "../particle/atom.txt"}}},
    {"Potential", {{"Model", "LennardJonesCutShift"}}},
    //{"Potential", {{"Model", "LennardJones"}}},
    //{"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "0.2"},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Log", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e2)}, {"tolerance", "1e-8"}}},
    {"ProfileCPU", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj_prof.csv"}}},
    {"Tune", {{}}},
    {"Run", {{"until_num_particles", "50"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Run", {{"num_trials", str(1e3)}}},
    {"Remove", {{"name", "Tune"}}},
    {"Run", {{"num_trials", str(1e3)}}},
    {"WriteCheckpoint", {{}}},
  }}, true);
  EXPECT_EQ(mc->random().class_name(), "RandomModulo");
  EXPECT_EQ(2, mc->configuration().num_particle_types());
  EXPECT_EQ(1, mc->system().unoptimized().num());
  //EXPECT_EQ(2, mc->system().unoptimized().num());
  EXPECT_EQ(1.2, mc->thermo_params().beta());
  EXPECT_EQ(1, mc->trials().num());
  EXPECT_EQ("TrialTranslate", mc->trial(0).class_name());
  //EXPECT_EQ("TrialAdd", mc->trial(1).class_name());
  EXPECT_EQ(1, mc->num_modifiers());
  EXPECT_EQ("CheckEnergy", mc->modify(0).class_name());
  EXPECT_LT(100, mc->trials().num_attempts());
  EXPECT_EQ(50, mc->configuration().num_particles());
  EXPECT_NE(0, mc->trials().num());
  mc->run(MakeRemove({{"all_trials", "true"}}));
  EXPECT_EQ(0, mc->trials().num());
  EXPECT_NE(0, mc->num_analyzers());
  mc->run(MakeRemove({{"all_analyzers", "true"}}));
  EXPECT_EQ(0, mc->num_analyzers());
  EXPECT_NE(0, mc->num_modifiers());
  mc->run(MakeRemove({{"all_modifiers", "true"}}));
  EXPECT_EQ(0, mc->num_modifiers());
}

TEST(MonteCarlo, gen_5_spce_in_triclinic) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"side_length", "22,22,22"},
      {"xy", "4"}, {"yz", "4"}, {"xz", "4"},
      {"particle_type0", "../particle/spce.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "5"}}},
  }}, true);
  FileXYZ().write_for_vmd("tmp/spce_triclinic.xyz", mc->configuration());
}

TEST(MonteCarlo, group_in_arglist) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type", "lj:../particle/lj.txt"}, {"group0", "first"}, {"first_particle_index", "0"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}, {"group", "first"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/ljrst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }}, true);
  EXPECT_EQ(2, mc->configuration().num_groups());
  EXPECT_NEAR(-2.060346185437E+00, mc->configuration().particle(2).site(0).position().coord(0), NEAR_ZERO);
  EXPECT_TRUE(std::abs(mc->configuration().particle(0).site(0).position().coord(0)-1.077169909511E+00)>1e-8);
}

TEST(MonteCarlo, two_configs) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "1234"}}},
    {"Configuration", {{"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type0", "../particle/lj.txt"}, {"group0", "first"}, {"first_particle_index", "0"}}},
    {"Configuration", {{"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"},
      {"particle_type0", "../particle/lj.txt"}, {"group0", "first"}, {"first_particle_index", "0"}}},
    {"Potential", {{"Model", "LennardJones"}, {"configuration_index", "0"}}},
    {"Potential", {{"Model", "LennardJones"}, {"configuration_index", "1"}}},
    {"ThermoParams", {{"beta", "100.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"configuration_index", "0"}}},
    {"TrialTranslate", {{"configuration_index", "1"}}},
    {"TrialTransfer", {{"particle_type", "0"}, {"configuration_index", "0"}}},
    {"TrialTransfer", {{"particle_type", "0"}, {"configuration_index", "1"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj0.xyz"}, {"configuration_index", "0"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj1.xyz"}, {"configuration_index", "1"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Run", {{"num_trials", "1e2"}}},
    {"ThermoParams", {{"beta", "100.2"}, {"chemical_potential", "-10."}}},
    {"Run", {{"num_trials", "1e2"}}},
  }}, true);
  EXPECT_EQ(2, mc->system().num_configurations());
  EXPECT_NE(mc->system().potential(0, 0).stored_energy(),
            mc->system().potential(0, 1).stored_energy());
//  INFO(mc->system().potential(0, 0).stored_energy() << " " <<
//       mc->system().potential(0, 1).stored_energy());
  EXPECT_NE(mc->system().potential(0, 1).stored_energy(), 0);
  EXPECT_GT(std::abs(mc->system().potential(0, 1).stored_energy() + 16.7903), 0.0001);
}

TEST(MonteCarlo, weight_per_number_fraction_no_particle_type) {
  TRY(
    auto mc = MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"},
                                                     {"particle_type1", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1,1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"tunable_param", "1."}}},
      {"TrialTransfer", {{"particle_type", "0"}}},
      {"TrialTransfer", {{"particle_type", "1"}}},
      {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.txt"}}},
      {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.xyz"}}},
      {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
      {"Tune", {{}}},
      {"Run", {{"num_trials", "1e2"}}},
    }}, true);
    CATCH_PHRASE("Trial::weight_per_number_fraction requires Trial::particle_type");
  );
}

TEST(MonteCarlo, weight_per_number_fraction_on_add) {
  TRY(
    auto mc = MakeMonteCarlo({{
      {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"},
                                                     {"particle_type1", "../particle/lj.txt"}}},
      {"Potential", {{"Model", "LennardJones"}}},
      {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1,1"}}},
      {"Metropolis", {{}}},
      {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "0"}}},
      {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}}},
      {"TrialTransfer", {{"weight_per_number_fraction", "1."}, {"particle_type", "0"}}},
      {"TrialTransfer", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}}},
      {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.txt"}}},
      {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.xyz"}}},
      {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
      {"Tune", {{}}},
      {"Run", {{"num_trials", "1e2"}}},
    }}, true);
    CATCH_PHRASE("weight_per_number_fraction is not implemented for PerturbAdd");
  );
}

TEST(MonteCarlo, weight_per_number_fraction) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"},
                                                   {"particle_type1", "../particle/lj.txt"},
                                                   {"particle_type2", "../particle/lj.txt"},
                                                   {"add_particles_of_type2", "1"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1,1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "0"}, {"number_fraction_exclude_type", "2"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}, {"number_fraction_exclude_type", "2"}}},
    {"TrialTransfer", {{"particle_type", "0"}}},
    {"TrialTransfer", {{"particle_type", "1"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }}, true);
  //INFO(feasst_str(mc->next_arg()));
  //EXPECT_EQ(mc->next_arg().first, "Run");
  for (int i = 0; i < 1e2; ++i) {
    mc->attempt();
    const int num0 = mc->configuration().num_particles_of_type(0);
    const int num1 = mc->configuration().num_particles_of_type(1);
    const int num2 = mc->configuration().num_particles_of_type(2);
    EXPECT_EQ(num2, 1);
    EXPECT_EQ(mc->trial(0).weight(), static_cast<double>(num0)/(num0+num1));
    EXPECT_EQ(mc->trial(1).weight(), static_cast<double>(num1)/(num0+num1));
    DEBUG("num0 " << num0 << " num1 " << num1 << " " << static_cast<double>(num0)/(num0+num1));
  }
}

}  // namespace feasst
