#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/hard_sphere.h"
#include "system/include/ideal_gas.h"
#include "system/include/dont_visit_model.h"
#include "system/include/visit_model_cell.h"
#include "system/include/lennard_jones.h"
#include "system/include/model_two_body_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/run.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/energy.h"
#include "steppers/include/num_particles.h"
#include "charge/include/trial_transfer_multiple.h"
#include "charge/include/check_net_charge.h"
#include "charge/include/slab_correction.h"
#include "charge/include/charge_screened.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"

namespace feasst {

/**
  Attempt to reproduce the average energy reported in

  https://doi.org/10.1063/1.476834

  Note that the Ewald vector cutoff definition may be slightly different.
 */
TEST(MonteCarlo, spce_nvt_VERY_LONG) {
  const int trials_per = 1e5;
  MonteCarlo mc;
  mc.set(spce({
    {"physical_constants", "CODATA2010"},
    {"cubic_side_length", "24.8586887"},
    {"alpha", str(5.6/24.8586887)},
    {"kmax_squared", "38"},
    {"xyz_file", "../plugin/charge/test/data/spce_sample_config_hummer_eq.xyz"}}));
  mc.set(MakeThermoParams({{"beta", str(1/kelvin2kJpermol(298, mc.configuration()))}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "0.2"}}));
  // without erfc table EXPECT_NEAR(mc.criteria().current_energy(), -24027.470339718111, 1e-10);
  EXPECT_NEAR(mc.criteria().current_energy(), -24027.470338455631, 1e-3);
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(trials_per)}, {"output_file", "tmp/spce_nvt"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(trials_per)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeTune());
  // mc.seek_num_particles(512);
  // FileXYZ().write("spce_sample_config_hummer_eq.xyz", mc.configuration());
  mc.add(MakeCheckProperties({{"trials_per_update", str(trials_per)}}));
  mc.add(MakeCPUTime({{"trials_per_write", str(5*trials_per)}}));
  // mc.attempt(1e6);
  auto energy = MakeEnergy({{"trials_per_write", str(trials_per)},
                            {"trials_per_update", "1"},
                            {"output_file", "tmp/spce_nvt_energy.txt"}});
  mc.add(energy);
  mc.attempt(1e6);
  test_serialize_unique(mc);
  INFO("energy: " << energy->energy().str());
  const double num = mc.configuration().num_particles();
  EXPECT_NEAR(-46.82*num,
              energy->energy().average(),
              8*std::sqrt(std::pow(energy->energy().block_stdev(), 2) +
                          std::pow(0.02*num, 2)));
  // FileXYZ().write("spce_sample_config_hummer.xyz", mc.configuration());
}

TEST(MonteCarlo, spce_gce_LONG) {
  const int trials_per = 1e2;
  MonteCarlo mc;
  mc.set(spce({{"alpha", str(5.6/24.8586887)}, {"kmax_squared", "38"}, {"cubic_side_length", str(24.8586887)}}));
  { const double sigma = mc.configuration().model_params().select("sigma").value(0);
    INFO("sigma " << sigma);
    mc.add_to_reference(MakePotential(
      MakeModelTwoBodyFactory(MakeLennardJones(),
                              MakeChargeScreened()),
      MakeVisitModelCell({{"min_length", str(sigma)}})));
  }
  const double beta = 1/kelvin2kJpermol(525);
  mc.set(MakeThermoParams({
    {"beta", str(beta)},
    {"chemical_potential", str(-8.14/beta)}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "0.2"}}));
  mc.add(MakeTrialTransfer({
    {"weight", "1."},
    {"particle_type", "0"},
    {"num_steps", "4"},
    {"reference_index", "0"}}));
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(trials_per)}, {"output_file", "tmp/spce_gce"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(trials_per)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeTune());
  mc.add(MakeCheckProperties({{"trials_per_update", str(trials_per)}}));
  auto num = MakeNumParticles({{"trials_per_write", str(trials_per)},
                               {"output_file", "tmp/spce_gce_num.txt"}});
  mc.add(num);
  mc.attempt(1e5);
  INFO(num->num_particles().str());
  EXPECT_NEAR(num->num_particles().average(), 9, 4);
}

// Fast test to be run with valgrind
TEST(MonteCarlo, spce) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.set(spce({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}}));
  mc.get_system()->add(MakePotential(MakeSlabCorrection({{"dimension", "0"}})));
  const double beta = 1/kelvin2kJpermol(525);
  mc.set(MakeThermoParams({
    {"beta", str(beta)},
    {"chemical_potential", str(-8.14/beta)}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "0.2"}}));
  mc.add(MakeTrialTransfer({{"weight", "4."}, {"particle_type", "0"}}));
  //mc.add(MakeTrialAdd({{"weight", "4."}, {"particle_type", "0"}}));
  //mc.run(MakeRun({{"until_num_particles", "2"}}));
  //mc.run(MakeRemove({{"name", "TrialAdd"}}));
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(5e2)}, {"output_file", "tmp/spce"}}));
  //mc.add(MakeCheckEnergyAndTune({{"trials_per", "1"}, {"tolerance", str(1e-6)}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(5e2)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeTune());
  mc.attempt(1e3);
  //EXPECT_GT(mc.configuration().num_particles(), 0);
}

TEST(MonteCarlo, spce_NVT_BENCHMARK_LONG) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.set(spce({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}}));
  const double beta = 1/kelvin2kJpermol(525);
  mc.set(MakeThermoParams({
    {"beta", str(beta)},
    {"chemical_potential", str(-8.14/beta)}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "0.2"}}));
  mc.add(MakeTrialTransfer({{"weight", "4."}, {"particle_type", "0"}}));
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(1e4)}, {"output_file", "tmp/spce"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(1e4)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeTune());
  mc.attempt(2e5);
}

TEST(MonteCarlo, rpm) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937());
  //mc.set(MakeRandomMT19937({{"seed", "123"}})); WARN("temporary");
  mc.set(rpm({{"alpha", str(5.6/20)}, {"kmax_squared", "38"}, {"cubic_side_length", "20"}}));
  //mc.set(rpm({{"cubic_side_length", "20"}, {"kmax_squared", "3"}})); WARN("temp");
  { Configuration * config = mc.get_system()->get_configuration();
    config->add_particle_of_type(0);
    config->add_particle_of_type(1);
    config->update_positions({{0., 0., 0.}, {1.01, 0., 0.}});
  }
  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
  const int trials_per = 1e0;
  mc.set(MakeThermoParams({
    {"beta", "0.02"},
    {"chemical_potential", "-509,-509"}}));
  mc.set(MakeMetropolis());
  EXPECT_NEAR(-0.99036730859815814, mc.criteria().current_energy(), 1e-9);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialTransferMultiple({
    {"weight", "1."},
    {"particle_types", "0,1"},
    {"reference_index", "0"},
  }));
  mc.add(MakeCheckProperties({{"trials_per_update", str(trials_per)}}));
  mc.add(MakeCheckPhysicality({{"trials_per_update", str(trials_per)}}));
  mc.add(MakeCPUTime({{"trials_per_write", str(5*trials_per)}}));
  mc.add(MakeCheckNetCharge());
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(trials_per)}, {"output_file", "tmp/rpm"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(trials_per)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeTune());
  mc.attempt(1e3);
  test_serialize_unique(mc);
}

TEST(MonteCarlo, spcearglist) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "20"},
                       {"particle_type0", "../particle/spce.txt"},
                       {"particle_type1", "../plugin/charge/particle/rpm_plus.txt"}}},
    {"Potential", {{"VisitModel", "Ewald"}, {"alpha", str(5.6/20)}, {"kmax_squared", "38"}}},
    {"Potential", {{"Model", "ModelTwoBodyFactory"}, {"model0", "LennardJones"}, {"model1", "ChargeScreened"}, {"VisitModel", "VisitModelCutoffOuter"}, {"erfc_table_size", "2e4"}}},
    {"Potential", {{"Model", "ChargeScreenedIntra"}, {"VisitModel", "VisitModelBond"}}},
    {"Potential", {{"Model", "ChargeSelf"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10,10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "0.2"},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Log", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e2)}, {"tolerance", "1e-8"}}},
    {"Tune", {{}}},
    {"Run", {{"until_num_particles", "50"}}},
    {"ThermoParams", {{"beta", "1.2"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Run", {{"num_trials", str(1e3)}}},
    {"Remove", {{"name", "Tune"}}},
    {"Run", {{"num_trials", str(1e3)}}},
    {"WriteCheckpoint", {{}}},
  }}, true);
  EXPECT_EQ("ModelTwoBodyFactory", mc->system().potential(1).model().class_name());
}

TEST(MonteCarlo, spce_npt) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "40"},
                       {"particle_type0", "../particle/spce.txt"}}},
    {"Potential", {{"VisitModel", "Ewald"}, {"alpha", str(5.6/20)}, {"kmax_squared", "38"}}},
    //{"Potential", {{"Model", "ModelTwoBodyFactory"}, {"model0", "LennardJones"}, {"model1", "ChargeScreened"}, {"erfc_table_size", "2e4"}}},
    {"Potential", {{"Model", "ModelTwoBodyFactory"}, {"model0", "LennardJones"}, {"model1", "ChargeScreened"}, {"VisitModel", "VisitModelCutoffOuter"}, {"erfc_table_size", "2e4"}}},
    {"Potential", {{"Model", "ChargeScreenedIntra"}, {"VisitModel", "VisitModelBond"}}},
    {"Potential", {{"Model", "ChargeSelf"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "0.01"}, {"chemical_potential", "10"}, {"pressure", "0.1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "0.2"},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialParticlePivot", {{"particle_type", "0"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    //{"Log", {{"trials_per_write", "1"}, {"output_file", "tmp/spce_npt.txt"}}},
    //{"Movie", {{"trials_per_write", "1"}, {"output_file", "tmp/spce_npt.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", "1e2"}, {"tolerance", "1e-8"}}},
    {"Tune", {{}}},
    {"Run", {{"until_num_particles", "10"}}},
    {"ThermoParams", {{"beta", "0.02"}, {"pressure", "0.01"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"TrialVolume", {{"tunable_target_acceptance", "0.5"}}},
    {"Run", {{"num_trials", str(1e3)}}},
    {"Remove", {{"name", "Tune"}}},
    {"Run", {{"num_trials", str(1e3)}}},
    {"WriteCheckpoint", {{}}},
  }}, true);
  //EXPECT_EQ("ModelTwoBodyTable", mc->system().potential(1).model().class_name());
}

}  // namespace feasst
