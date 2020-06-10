#include "utils/test/utils.h"
#include "system/include/utils.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "growth_expanded/include/trial_growth_expanded.h"
#include "growth_expanded/include/macrostate_growth_expanded.h"

namespace feasst {

TEST(MonteCarlo, TrialGrowthExpanded) {
  const std::string data = "forcefield/data.dimer";
  MonteCarlo mc;
  mc.set(lennard_jones({{"particle", data}}));
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e2)}, {"file_name", "tmp/growth"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e2)}}));
  mc.set(MakeMetropolis({{"beta", "1"}, {"chemical_potential", "1."}}));
  SeekNumParticles(4).with_trial_add().run(&mc);
  mc.set(MakeFlatHistogram(
    MakeMacrostateGrowthExpanded(Histogram({{"width", "0.5"}, {"max", "10"}})),
    MakeTransitionMatrix({{"min_sweeps", "10"}}),
    // MakeWangLandau({{"min_flatness", "1"}}),
    { {"beta", str(1./1.5)},
      {"chemical_potential", "-6.952321"}}));

  mc.add(MakeTrialGrowthExpanded(build_(1, data), build_(2, data)));
  mc.add(MakeCriteriaWriter({{"steps_per_write", str(int(1e6))}}));
  mc.add(MakeCriteriaUpdater({{"steps_per_write", str(int(1e6))}}));
  EXPECT_EQ(2, mc.trials().num());
//  for (int i = 0; i < 50; ++i) {

// new bug introduced where,when adding new particle, sites need to be set as
// non physical for growth expanded
//  mc.attempt(1e2);

//    INFO(mc.criteria()->write());
//  }
  test_serialize(mc);
//  std::stringstream ss;
//  mc.serialize(ss);
//  INFO(ss.str());
}

}  // namespace feasst
