#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "steppers/include/movie.h"
#include "steppers/include/log.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/energy.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "steppers/include/criteria_writer.h"
#include "ewald/test/system_example.h"
#include "ewald/include/check_net_charge.h"
#include "egce/include/a_equal_or_one_more_than_b.h"

namespace feasst {

TEST(MonteCarlo, rpm_egce_fh) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  {
    Configuration config(MakeDomain({{"cubic_box_length", "12"}}), {
      {"particle_type0", "../plugin/ewald/forcefield/data.rpm_plus"},
      {"particle_type1", "../plugin/ewald/forcefield/data.rpm_minus"}
    });
    const double rcut = 4.891304347826090;
    config.set_model_param("cutoff", 0, rcut);
    config.set_model_param("cutoff", 1, rcut);
    config.set_model_param("charge", 0, 1./std::sqrt(CODATA2018().charge_conversion()));
    config.set_model_param("charge", 1, -1./std::sqrt(CODATA2018().charge_conversion()));
    mc.add(config);
  }
  mc.add(Potential(
    MakeEwald({{"kmax_squared", "38"},
               {"alpha", str(6.87098396396261/mc.configuration().domain()->min_side_length())}})));
  mc.add(Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                            MakeChargeScreened()})));
  mc.add(Potential(MakeChargeSelf()));
  const double temperature = 0.047899460618081;
  const double beta_mu = -13.94;
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "20"}, {"min", "0"}}),
      MakeAEqualOrOneMoreThanB()),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    { {"beta", str(1/temperature)},
      {"chemical_potential0", str(beta_mu*temperature)},
      {"chemical_potential1", str(beta_mu*temperature)}}
  );
  mc.set(criteria);
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.1"}}));
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "0"}});
  add_trial_transfer(&mc, {{"weight", "1."}, {"particle_type", "1"}});
  const int steps_per = 1e6;
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/rpm_egce_crit.txt"},
  }));
  mc.add(MakeMovie({{"file_name", "tmp/rpm_egce.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/rpm_egce_log.txt"}, {"steps_per", str(steps_per)}}));
  // mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(1e6)}}));
  // mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeEnergy({
    {"file_name", "wlmc_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"},
  }));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckNetCharge({{"maximum", "1."}, {"minimum", str(-NEAR_ZERO)}}));
  mc.attempt(1e3);
  // mc.run_until_complete();
  test_serialize(mc);
}

}  // namespace feasst
