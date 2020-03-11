#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "system/include/dont_visit_model.h"
#include "system/include/model_two_body_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "flat_histogram/test/flat_histogram_test.h"
#include "math/include/histogram.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/check_physicality.h"
#include "ewald/include/ewald.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/charge_self.h"
#include "ewald/include/check_net_charge.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/trial_remove_multiple.h"

namespace feasst {

const double energy_av(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

TEST(FlatHistogram, order) {
  auto criteria = MakeFlatHistogram({
    {"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}
  });
  TRY(
    auto criteria2 = criteria;
    criteria2->set(MakeWangLandau({{"min_flatness", "20"}}));
    CATCH_PHRASE("set macrostate before bias");
  );
}

TEST(TransitionMatrix, args) {
  TRY(
    auto criteria = crit_fh(0);
    criteria->set(MakeTransitionMatrix());
    CATCH_PHRASE("key(min_sweeps) is required");
  );
}

TEST(WangLandau, args) {
  TRY(
    auto criteria = crit_fh(1);
    criteria->set(MakeWangLandau());
    CATCH_PHRASE("key(min_flatness) is required");
  );
}

TEST(MonteCarlo, FHMC) {
  // for (int crit_type = 0; crit_type < 1; ++crit_type) {
  // for (int crit_type = 1; crit_type < 2; ++crit_type) {
  for (int crit_type = 0; crit_type < 3; ++crit_type) {
    MonteCarlo mc;
    // mc.set(MakeRandomMT19937({{"seed", "default"}}));
    //mc.set(MakeRandomMT19937({{"seed", "1583434540"}}));
    mc_lj(&mc);
    // mc.seek_num_particles(4);
    add_trial_transfer(&mc, {{"particle_type", "0"}, {"weight", "0.25"}});
    auto crit = crit_fh(crit_type);
    INFO(crit->bias()->class_name());
    mc.set(crit);//crit_fh(crit_type));
    mc.add(MakeMovie({
      {"file_name", "tmp/wlmc_movie"},
      {"steps_per", str(1e4)},
      {"multistate", "true"},
    }));
    mc.add(MakeCriteriaUpdater({{"steps_per", str(1)}}));
    mc.add(MakeCriteriaWriter({
      {"steps_per", str(1e4)},
      {"file_name", "tmp/ljcrit.txt"},
    }));
    mc.add(MakeEnergy({
      {"file_name", "wlmc_energy"},
      {"steps_per_update", "1"},
      {"steps_per_write", str(1e4)},
      {"multistate", "true"},
    }));
    // mc.attempt(1e4);
    // mc.attempt(1e6); // note more than 1e4 steps required for TM
    mc.run_until_complete();
    // INFO(mc.criteria()->write());

    //MonteCarlo mc2 = test_serialize_no_comp(mc);
    test_serialize(mc);

    // compare with known values of lnpi
    const LnProbability * lnpi = &crit->bias()->ln_prob();
    EXPECT_NEAR(lnpi->value(0), -18.707570324988800000, 0.55);
    EXPECT_NEAR(lnpi->value(1), -14.037373358321800000, 0.55);
    EXPECT_NEAR(lnpi->value(2), -10.050312091655200000, 0.55);
    EXPECT_NEAR(lnpi->value(3), -6.458920624988570000, 0.55);
    EXPECT_NEAR(lnpi->value(4), -3.145637424988510000, 0.55);
    EXPECT_NEAR(lnpi->value(5), -0.045677458321876000, 0.55);

    // compare with known values of energy
    EXPECT_NEAR(energy_av(0, mc), 0, 1e-14);
    EXPECT_NEAR(energy_av(1, mc), -0.000605740233333333, 1e-8);
    EXPECT_NEAR(energy_av(2, mc), -0.030574223333333334, 0.02);
    EXPECT_NEAR(energy_av(3, mc), -0.089928316, 0.05);
    EXPECT_NEAR(energy_av(4, mc), -0.1784570533333333, 0.06);
    EXPECT_NEAR(energy_av(5, mc), -0.29619201333333334, 0.1);
    EXPECT_LE(mc.system().configuration().num_particles(), 5);
  }
}

TEST(MonteCarlo, rpm_fh) {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  { Configuration config(
      MakeDomain({{"cubic_box_length", "12"}}),
      {{"particle_type0", "../plugin/ewald/forcefield/data.rpm_plus"},
       {"particle_type1", "../plugin/ewald/forcefield/data.rpm_minus"}}
    );
    const double rcut = 4.891304347826090;
    config.set_model_param("cutoff", 0, rcut);
    config.set_model_param("cutoff", 1, rcut);
    config.set_model_param("charge", 0, 1./std::sqrt(CODATA2018().charge_conversion()));
    config.set_model_param("charge", 1, -1./std::sqrt(CODATA2018().charge_conversion()));
//    config.add_particle_of_type(0);
//    config.update_positions({{1.5, 0., 0.}});
//    config.add_particle_of_type(1);
//    config.update_positions({{0., 0., 0.}, {1.5, 0., 0.}});
    mc.add(config);
  }
  mc.add(Potential(
    MakeEwald({{"kmax_squared", "38"},
               {"alpha", str(6.87098396396261/mc.configuration().domain()->min_side_length())}})));
  mc.add(Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                            MakeChargeScreened()})));
  mc.add(Potential(MakeChargeSelf()));
  mc.add_to_reference(Potential(MakeDontVisitModel()));
  //const double temperature = 0.047899460618081;
  INFO("charge conversion " << CODATA2018().charge_conversion());
  const double temperature = 0.047899460618081;
  // const double temperature = 0.047899460618081*CODATA2018().charge_conversion();
  const double beta_mu = -13.94;
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      //Histogram({{"width", "1"}, {"max", "1"}, {"min", "0"}}),
      //Histogram({{"width", "1"}, {"max", "2"}, {"min", "1"}}),
      Histogram({{"width", "1"}, {"max", "2"}, {"min", "0"}}),
      {{"particle_type", "0"}}
    ),
    // MakeWangLandau({{"min_flatness", "15"}}),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    {{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)},
    }
  );
  mc.set(criteria);
  INFO("beta_mu " << criteria->beta_mu(0));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.1"}}));
  const argtype transfer_args = {
    {"weight", "1."},
    {"particle_type0", "0"},
    {"particle_type1", "1"},
    {"reference_index", "0"},
  };
  mc.add(MakeTrialAddMultiple(transfer_args));
  mc.add(MakeTrialRemoveMultiple(transfer_args));

  const int steps_per = 1e6;
  if (criteria->bias()->class_name() == "TransitionMatrix") {
    mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  }
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/rpmcrit.txt"},
  }));
  mc.add(MakeMovie({{"file_name", "tmp/rpm.xyz"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeLog({{"file_name", "tmp/rpm_log.txt"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  //mc.add(MakeTuner({{"steps_per", str(1e2)}}));
  mc.add(MakeEnergy({
    {"file_name", "wlmc_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"},
  }));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckPhysicality({{"steps_per", str(steps_per)}}));
  mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckEnergy({{"tolerance", "1e-8"}, {"steps_per", str(steps_per)}}));
  mc.add(MakeCheckNetCharge());
  mc.attempt(1e4);
  // mc.run_until_complete();

  test_serialize(mc);

  const LnProbability * lnpi = &criteria->bias()->ln_prob();
  for (int bin = 0; bin < lnpi->size(); ++bin) {
    INFO(bin << " " << lnpi->value(bin));
  }
}

}  // namespace feasst
