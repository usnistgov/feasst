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
  // for (int crit_type = 0; crit_type < 1; ++crit_type) {
  for (int crit_type = 0; crit_type < 2; ++crit_type) {
    MonteCarlo mc;
    // mc.set(MakeRandomMT19937({{"seed", "default"}}));
    mc_lj(&mc);
    // mc.seek_num_particles(4);
    add_trial_transfer(&mc, {{"particle_type", "0"}, {"weight", "0.25"}});
    mc.set(crit_fh(crit_type));
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

    INFO(mc2.analyzers().back()->analyzers().back()->accumulator().average());

    // compare with known values of lnpi and energy
    //const LnProbability * lnpi = &criteria->bias()->ln_macro_prob();
    //INFO(lnpi->value(0));
    INFO(energy_av(0, mc));
    EXPECT_LE(mc.system().configuration().num_particles(), 5);
  }
}

}  // namespace feasst
