#include "utils/test/utils.h"
#include "math/include/histogram.h"
#include "math/include/random_mt19937.h"
#include "system/include/potential.h"
#include "system/include/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/metropolis.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/wltm.h"
#include "beta_expanded/include/trial_beta.h"
#include "beta_expanded/include/macrostate_beta.h"

namespace feasst {

TEST(MonteCarlo, beta_expanded) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.set(lennard_jones());
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  const double beta_min = 0.8;
  const double beta_max = 1.2;
  const int beta_num = 5;
  const std::string delta_beta = str((beta_max - beta_min)/(beta_num-1));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."},{"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "10"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.set(MakeFlatHistogram({{"Macrostate", "MacrostateBeta"}, {"width", delta_beta}, {"max", str(beta_max)}, {"min", str(beta_min)},
    {"Bias", "WLTM"}, {"collect_flatness", "18"}, {"min_flatness", "22"}, {"min_sweeps", "10"}}));
  mc.add(MakeTrialBeta({{"fixed_beta_change", delta_beta}}));
  const std::string steps_per(str(1e4));
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/lj_beta"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", steps_per}}));
  mc.add(MakeCriteriaUpdater({{"steps_per", steps_per}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", steps_per},
    {"file_name", "tmp/lj_beta_crit.txt"},
    {"file_name_append_phase", "true"}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/lj_beta_energy"},
    {"file_name_append_phase", "true"},
    {"steps_per_update", "1"},
    {"steps_per_write", steps_per},
    {"multistate", "true"}}));
  mc.attempt(5e4);
}

}  // namespace feasst
