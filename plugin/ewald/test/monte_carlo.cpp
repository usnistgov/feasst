#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "system/include/ideal_gas.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/cpu_time.h"
#include "ewald/test/system_example.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/check_net_charge.h"

namespace feasst {

TEST(MonteCarlo, spce) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.add(Configuration(
    MakeDomain({{"cubic_box_length", "24.8586887"}}),
    {{"particle_type", "../forcefield/data.spce"}}
  ));
  mc.add(Potential(
    MakeEwald({{"kmax_squared", "3"},
               {"alpha", str(5.6/mc.configuration().domain()->min_side_length())}})
  ));
  mc.add(Potential(MakeModelTwoBodyFactory({MakeLennardJones(),
                                            MakeChargeScreened()})));
  mc.add(Potential(MakeChargeScreenedIntra(), MakeVisitModelBond()));
  mc.add(Potential(MakeChargeSelf()));
  mc.add(Potential(MakeLongRangeCorrections()));
  const int steps_per = 1e2;
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-9"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "1."}}));
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  mc.add(MakeMovie({{"file_name", "tmp/spce.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/spce_log.txt"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(1e2)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.attempt(1e3);
  test_serialize(mc);
}

TEST(MonteCarlo, rpm) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  { Configuration config(
      MakeDomain({{"cubic_box_length", "20"}}),
      {{"particle_type0", "../plugin/ewald/forcefield/data.rpm_plus"},
       {"particle_type1", "../plugin/ewald/forcefield/data.rpm_minus"}}
    );
    config.add_particle_of_type(0);
    config.add_particle_of_type(1);
    config.update_positions({{0., 0., 0.}, {1.01, 0., 0.}});
    mc.add(config);
  }
  mc.add(Potential(
    MakeEwald({{"kmax_squared", "38"},
               {"alpha", str(5.6/mc.configuration().domain()->min_side_length())}})));
  mc.add(Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                            MakeChargeScreened()})));
  mc.add(Potential(MakeChargeSelf()));
  mc.add_to_reference(Potential(MakeDontVisitModel()));
  const int steps_per = 1e2;
  mc.set(MakeMetropolis({
    {"beta", "0.02"},
    {"chemical_potential0", "-509"},
    {"chemical_potential1", "-509"},
  }));
  EXPECT_NEAR(-0.99036730859815814,
    mc.criteria()->current_energy()/CODATA2018().charge_conversion(),
    1e-9);
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  const argtype transfer_args = {
    {"weight", "1."},
    {"particle_type0", "0"},
    {"particle_type1", "1"},
    {"reference_index", "0"},
  };
  mc.add(MakeTrialAddMultiple(transfer_args));
  mc.add(MakeTrialRemoveMultiple(transfer_args));
  mc.add(MakeMovie({{"file_name", "tmp/rpm.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/rpm_log.txt"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(1e2)}}));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckPhysicality({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeCheckNetCharge());
  mc.attempt(1e3);
  test_serialize(mc);
}

}  // namespace feasst
