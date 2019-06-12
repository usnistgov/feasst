#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/bias_wang_landau.h"
#include "flat_histogram/include/bias_transition_matrix.h"
#include "math/include/histogram.h"
#include "steppers/include/energy.h"

namespace feasst {

TEST(MonteCarlo, WLMC) {
  for (int crit_type = 0; crit_type < 2; ++crit_type) {
    MonteCarlo mc = mc_lj();
    mc.add(MakeTrialAdd({{"weight", "0.25"}}));
    mc.add(MakeTrialRemove({{"weight", "0.25"}}));
    { auto criteria = MakeCriteriaFlatHistogram({
        {"beta", str(1./1.5)},
        {"chemical_potential", "-2.352321"}
      });
      try {
        auto criteria2 = criteria;
        criteria2->set(MakeBiasWangLandau({{"min_flatness", "20"}}));
        CATCH_PHRASE("set macrostate before bias");
      }
      criteria->set(MakeMacrostateNumParticles(Histogram({
        {"width", "1"},
        {"max", "5"},
      })));
      if (crit_type == 0) {
        try {
          auto criteria2 = criteria;
          criteria2->set(MakeBiasTransitionMatrix());
          CATCH_PHRASE("key(min_sweeps) is required");
        }
        criteria->set(MakeBiasTransitionMatrix({
          {"min_sweeps", "20"},
          {"num_steps_to_update", "1"},
        }));
      } else {
        try {
          auto criteria2 = criteria;
          criteria2->set(MakeBiasWangLandau());
          CATCH_PHRASE("key(min_flatness) is required");
        }
        criteria->set(MakeBiasWangLandau({{"min_flatness", "20"}}));
      }
      mc.set(criteria);
    }
    mc.add(MakeMovie({
      {"file_name", "tmp/wlmc_movie"},
      {"steps_per", "1000"},
      {"multistate", "true"},
    }));
    mc.add(MakeEnergy({
      {"file_name", "tmp/wlmc_energy"},
      {"steps_per_update", "1"},
      {"steps_per_write", "1"},
      {"multistate", "true"},
    }));
    mc.attempt(1e4);
    // mc.attempt(1e6); // note more than 1e4 steps required for TM
    // mc.run_until_complete();
    INFO(mc.criteria()->write());

    MonteCarlo mc2 = test_serialize_no_comp(mc);
  }
}

}  // namespace feasst
