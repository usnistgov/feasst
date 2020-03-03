#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "flat_histogram/test/flat_histogram_test.h"
#include "math/include/histogram.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"

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
  for (int crit_type = 0; crit_type < 1; ++crit_type) {
  // for (int crit_type = 0; crit_type < 2; ++crit_type) {
    MonteCarlo mc;
    // mc.set(MakeRandomMT19937({{"seed", "default"}}));
    mc_lj(&mc);
    // mc.seek_num_particles(4);
    add_trial_transfer(&mc, {{"particle_type", "0"}, {"weight", "0.25"}});
    auto crit = crit_fh(crit_type);
    mc.set(crit);//crit_fh(crit_type));
    mc.add(MakeMovie({
      {"file_name", "tmp/wlmc_movie"},
      {"steps_per", str(1e4)},
      {"multistate", "true"},
    }));
    if (crit_type == 0) {
      mc.add(MakeCriteriaUpdater({{"steps_per", str(1)}}));
    }
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
    INFO(mc.criteria()->write());

    MonteCarlo mc2 = test_serialize_no_comp(mc);

    // compare with known values of lnpi
    const LnProbability * lnpi = &crit->bias()->ln_prob();
    EXPECT_NEAR(lnpi->value(0), -18.707570324988800000, 0.4);
    EXPECT_NEAR(lnpi->value(1), -14.037373358321800000, 0.4);
    EXPECT_NEAR(lnpi->value(2), -10.050312091655200000, 0.4);
    EXPECT_NEAR(lnpi->value(3), -6.458920624988570000, 0.4);
    EXPECT_NEAR(lnpi->value(4), -3.145637424988510000, 0.4);
    EXPECT_NEAR(lnpi->value(5), -0.045677458321876000, 0.4);

    // compare with known values of energy
    EXPECT_NEAR(energy_av(0, mc), 0, 1e-14);
    EXPECT_NEAR(energy_av(1, mc), -0.000605740233333333, 1e-8);
    EXPECT_NEAR(energy_av(2, mc), -0.030574223333333334, 0.02);
    EXPECT_NEAR(energy_av(3, mc), -0.089928316, 0.03);
    EXPECT_NEAR(energy_av(4, mc), -0.1784570533333333, 0.04);
    EXPECT_NEAR(energy_av(5, mc), -0.29619201333333334, 0.1);
    EXPECT_LE(mc.system().configuration().num_particles(), 5);
  }
}

}  // namespace feasst
