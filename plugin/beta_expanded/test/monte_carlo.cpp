#include "utils/test/utils.h"
#include "math/include/histogram.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "system/include/thermo_params.h"
#include "system/include/potential.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/remove_trial.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/metropolis.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/wltm.h"
#include "beta_expanded/include/trial_beta.h"
#include "beta_expanded/include/macrostate_beta.h"

namespace feasst {

TEST(MonteCarlo, beta_expanded) {
  const double beta_min = 0.8;
  const double beta_max = 1.2;
  const int beta_num = 5;
  const std::string delta_beta = str((beta_max - beta_min)/(beta_num-1));
  const std::string trials_per(str(1e4));
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."},{"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "10"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateBeta"}, {"width", delta_beta}, {"max", str(beta_max)}, {"min", str(beta_min)},
      {"Bias", "WLTM"}, {"collect_flatness", "18"}, {"min_flatness", "22"}, {"min_sweeps", "10"}}},
    {"TrialBeta", {{"fixed_beta_change", delta_beta}}},
//    {"LogAndMovie", {{"trials_per_write", trials_per}, {"output_file", "tmp/lj_beta"}}},
    {"CheckEnergy", {{"trials_per_update", trials_per}}},
    {"Tune", {{}}},
    {"CriteriaUpdater", {{"trials_per_update", trials_per}}},
    {"CriteriaWriter", {{"trials_per_write", trials_per},
      {"output_file", "tmp/lj_beta_crit.txt"},
      {"output_file_append_phase", "true"}}},
    {"Energy", {
      {"output_file", "tmp/lj_beta_energy"},
      {"output_file_append_phase", "true"},
      {"trials_per_update", "1"},
      {"trials_per_write", trials_per},
      {"multistate", "true"}}}
  }}, true);
  mc->attempt(5e4);
}

//MonteCarlo sweeptest(const int min_sweeps, const int beta_num) {
//  MonteCarlo mc;
//  mc.add(MakeConfiguration({{"cubic_side_length", "8"}, {"particle_type0", "../particle/lj.txt"}}));
//  mc.add(MakePotential(MakeDontVisitModel()));
//  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
//  const double beta_min = 0.8;
//  const double beta_max = 1.2;
////  const int beta_num = 50;
//  const std::string delta_beta = str((beta_max - beta_min)/(beta_num-1));
//  //INFO("delta_beta " << delta_beta);
//  mc.set(MakeFlatHistogram({{"Macrostate", "MacrostateBeta"}, {"width", delta_beta}, {"max", str(beta_max)}, {"min", str(beta_min)},
//    {"Bias", "TransitionMatrix"}, {"min_sweeps", str(min_sweeps)}, {"min_visits", "1"}, {"average_visits", "1000"}, {"new_sweep", "1"}}));
//  mc.add(MakeTrialBeta({{"fixed_beta_change", delta_beta}}));
//  const std::string trials_per(str(1e6));
//  mc.add(MakeLogAndMovie({{"trials_per", trials_per}, {"output_file", "tmp/lj_beta"}}));
//  mc.add(MakeCheckEnergyAndTune({{"trials_per", trials_per}}));
//  mc.add(MakeCriteriaUpdater({{"trials_per", "1"}}));
//  //mc.add(MakeCriteriaUpdater({{"trials_per", trials_per}}));
//  mc.add(MakeCriteriaWriter({
//    {"trials_per", trials_per},
//    {"output_file", "tmp/lj_beta_crit.txt"},
//    {"output_file_append_phase", "true"}}));
//  mc.add(MakeEnergy({
//    {"output_file", "tmp/lj_beta_energy"},
//    {"output_file_append_phase", "true"},
//    {"trials_per_update", "1"},
//    {"trials_per_write", trials_per},
//    {"multistate", "true"}}));
//  return mc;
//}
//
//TEST(MonteCarlo, sweeping_window_size_dependence_LONG) {
//  for (int num : {2, 5, 10, 20, 50, 100, 200, 500}) {
//    Accumulator num_trials;
//    for (int i = 0; i < 3; ++i) {
//      MonteCarlo mc = sweeptest(1e4, num);
//      mc.run_until_complete();
//      num_trials.accumulate(mc.trials().num_attempts());
//    }
//    INFO(num << " " << num_trials.average() << " " << num_trials.std() << " t/n " << num_trials.average()/num);
//  }
//}

}  // namespace feasst
