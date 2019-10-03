
#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "flat_histogram/include/bias_transition_matrix.h"
#include "flat_histogram/include/bias_wang_landau.h"
#include "steppers/include/criteria_writer.h"
#include "growth_expanded/include/trial_growth_expanded.h"
#include "growth_expanded/include/macrostate_growth_expanded.h"

namespace feasst {

TEST(MonteCarlo, TrialGrowthExpanded) {
  seed_random_by_date();
  // seed_random(1563999738);
  const std::string data = "../forcefield/data.dimer";
  MonteCarlo mc;
  mc_lj(&mc, 8, data, 1e2, true);
  mc.set(MakeCriteriaMetropolis({{"beta", "1"}, {"chemical_potential", "1."}}));
  mc.seek_num_particles(4);

  { auto criteria = MakeCriteriaFlatHistogram({{"beta", str(1./1.5)},
      {"chemical_potential", "-6.952321"}});
    criteria->set(MakeMacrostateGrowthExpanded(
      Histogram({{"width", "0.5"}, {"max", "10"}})//,
//      {{"soft_max", "10"}}
    ));
    criteria->set(MakeBiasTransitionMatrix({
      {"min_sweeps", "10"},
      {"num_steps_to_update", str(1e6)}, // fast
    }));
    // criteria->set(MakeBiasWangLandau({{"min_flatness", "1"}}));
    mc.set(criteria);
  }

  mc.add(MakeTrialGrowthExpanded(build_(1, data), build_(2, data)));
  mc.add(MakeCriteriaWriter({{"steps_per_write", str(int(1e6))}}));
  EXPECT_EQ(2, mc.trials().num_trials());
//  for (int i = 0; i < 50; ++i) {
  mc.attempt(1e2);
//    INFO(mc.criteria()->write());
//  }
  test_serialize(mc);
//  std::stringstream ss;
//  mc.serialize(ss);
//  INFO(ss.str());
}

}  // namespace feasst
