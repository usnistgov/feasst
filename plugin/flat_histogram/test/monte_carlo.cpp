#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/bias_wang_landau.h"
#include "flat_histogram/include/bias_transition_matrix.h"
#include "math/include/histogram.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"

namespace feasst {

const double energy_av(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

TEST(CriteriaFlatHistogram, order) {
  auto criteria = MakeCriteriaFlatHistogram({
    {"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}
  });
  try {
    auto criteria2 = criteria;
    criteria2->set(MakeBiasWangLandau({{"min_flatness", "20"}}));
    CATCH_PHRASE("set macrostate before bias");
  }
}

std::shared_ptr<CriteriaFlatHistogram> crit_fh(const int crit_type) {
  auto criteria = MakeCriteriaFlatHistogram({
    {"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}
  });
  criteria->set(MakeMacrostateNumParticles(
    Histogram({{"width", "1"}, {"max", "5"}}),
    {{"soft_max", "5"}}
  ));
  if (crit_type == 0) {
    criteria->set(MakeBiasTransitionMatrix({
      {"min_sweeps", "10"},
      // {"num_steps_to_update", str(1e5)},  // benchmark 1.7 seconds
      {"num_steps_to_update", "1"}, // fast
    }));
  } else {
    criteria->set(MakeBiasWangLandau({{"min_flatness", "1"}}));
  }
  return criteria;
}

TEST(BiasTransitionMatrix, args) {
  try {
    auto criteria = crit_fh(0);
    criteria->set(MakeBiasTransitionMatrix());
    CATCH_PHRASE("key(min_sweeps) is required");
  }
}

TEST(BiasWangLandau, args) {
  try {
    auto criteria = crit_fh(1);
    criteria->set(MakeBiasWangLandau());
    CATCH_PHRASE("key(min_flatness) is required");
  }
}

TEST(MonteCarlo, FHMC) {
  seed_random();
  // seed_random_by_date();
  // seed_random(1560889975);
  for (int crit_type = 0; crit_type < 1; ++crit_type) {
  // for (int crit_type = 0; crit_type < 2; ++crit_type) {
    MonteCarlo mc;
    mc_lj(&mc);
    mc.add(MakeTrialAdd({{"particle_type", "0"}, {"weight", "0.25"}}));
    mc.add(MakeTrialRemove({{"weight", "0.25"}}));
    mc.set(crit_fh(crit_type));
    mc.add(MakeMovie({
      {"file_name", "tmp/wlmc_movie"},
      {"steps_per", str(1e4)},
      {"multistate", "true"},
    }));
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
    //const LnProbabilityDistribution * lnpi = &criteria->bias()->ln_macro_prob();
    //INFO(lnpi->value(0));
    INFO(energy_av(0, mc));
    EXPECT_LE(mc.system().configuration().num_particles(), 5);
  }
}

}  // namespace feasst
