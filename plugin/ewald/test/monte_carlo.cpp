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
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/utils.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/utils.h"
#include "steppers/include/energy.h"
#include "steppers/include/num_particles.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/check_net_charge.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/utils.h"

namespace feasst {

/**
  Attempt to reproduce the average energy reported in

  https://doi.org/10.1063/1.476834

  Note that the Ewald vector cutoff definition may be slightly different.
 */
TEST(MonteCarlo, spce_nvt_LONG) {
  const int steps_per = 1e5;
  MonteCarlo mc;
  mc.set(spce({
    {"physical_constants", "CODATA2010"},
    {"cubic_box_length", "24.8586887"},
    {"alphaL", "5.6"},
    {"kmax_squared", "38"},
    {"xyz_file", install_dir() + "/plugin/ewald/test/data/spce_sample_config_hummer_eq.xyz"}}));
  mc.set(MakeMetropolis({
    {"beta", str(1/kelvin2kJpermol(298, mc.configuration()))}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
  EXPECT_NEAR(mc.criteria().current_energy(), -24027.470339718111, 1e-10);
  // add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  add_common_steppers(&mc, {{"steps_per", str(steps_per)},
                            {"file_append", "tmp/spce_nvt"},
                            {"tolerance", str(1e-6)}});
  // mc.seek_num_particles(512);
  // FileXYZ().write("spce_sample_config_hummer_eq.xyz", mc.configuration());
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  // mc.attempt(1e6);
  auto energy = MakeEnergy({{"steps_per_write", str(steps_per)},
                            {"steps_per_update", "1"},
                            {"file_name", "tmp/spce_nvt_energy.txt"}});
  mc.add(energy);
  mc.attempt(1e6);
  test_serialize(mc);
  INFO("energy: " << energy->energy().str());
  const double num = mc.configuration().num_particles();
  EXPECT_NEAR(-46.82*num,
              energy->energy().average(),
              8*std::sqrt(std::pow(energy->energy().block_stdev(), 2) +
                          std::pow(0.02*num, 2)));
  // FileXYZ().write("spce_sample_config_hummer.xyz", mc.configuration());
}

TEST(MonteCarlo, spce_gce_LONG) {
  const int steps_per = 1e2;
  MonteCarlo mc;
  mc.set(spce());
  { const double sigma = mc.configuration().model_params().sigma().value(0);
    INFO("sigma " << sigma);
    mc.get_system()->get_configuration()->get_domain()->init_cells(sigma);
    mc.add_to_reference(Potential(
      MakeModelTwoBodyFactory({MakeLennardJones(),
                               MakeChargeScreened()}),
      MakeVisitModelCell()));
  }
  mc.set(MakeMetropolis({
    {"beta", str(1/kelvin2kJpermol(525))},
    {"chemical_potential", "-35.294567543492"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
  add_trial_transfer(&mc,
    { {"weight", "1."},
      {"particle_type", "0"},
      {"num_steps", "4"},
      {"reference_index", "0"}});
  add_common_steppers(&mc, {{"steps_per", str(steps_per)},
                            {"file_append", "tmp/spce_gce"},
                            {"tolerance", str(1e-6)}});
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  auto num = MakeNumParticles({{"steps_per_write", str(steps_per)},
                               {"file_name", "tmp/spce_gce_num.txt"},
                               {"num_block", str(steps_per)}});
  mc.add(num);
  mc.attempt(1e5);
  INFO(num->num_particles().str());
  EXPECT_NEAR(num->num_particles().average(), 9, 2);
}

// Fast test to be run with valgrind
TEST(MonteCarlo, spce) {
  MonteCarlo mc;
  mc.set(spce());
  mc.set(MakeMetropolis({
    {"beta", str(1/kelvin2kJpermol(525))},
    {"chemical_potential", "-35.294567543492"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  add_common_steppers(&mc, {{"steps_per", str(5e2)},
                            {"file_append", "tmp/spce"},
                            {"tolerance", str(1e-6)}});
  mc.attempt(1e3);
}

TEST(MonteCarlo, rpm) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(rpm({{"cubic_box_length", "20"}}));
  { Configuration * config = mc.get_system()->get_configuration();
    config->add_particle_of_type(0);
    config->add_particle_of_type(1);
    config->update_positions({{0., 0., 0.}, {1.01, 0., 0.}});
  }
  mc.add_to_reference(Potential(MakeDontVisitModel()));
  const int steps_per = 1e2;
  mc.set(MakeMetropolis({
    {"beta", "0.02"},
    {"chemical_potential0", "-509"},
    {"chemical_potential1", "-509"},
  }));
  EXPECT_NEAR(-0.99036730859815814, mc.criteria().current_energy(), 1e-9);
  const argtype transfer_args = {
    {"weight", "1."},
    {"particle_type0", "0"},
    {"particle_type1", "1"},
    {"reference_index", "0"},
  };
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialAddMultiple(transfer_args));
  mc.add(MakeTrialRemoveMultiple(transfer_args));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckPhysicality({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckNetCharge());
  add_common_steppers(&mc, {{"steps_per", str(steps_per)},
                            {"file_append", "tmp/rpm"}});
  mc.attempt(1e3);
  test_serialize(mc);
}


}  // namespace feasst
