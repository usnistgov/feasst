
#ifndef FEASST_CORE_UI_BRIEF_H_
#define FEASST_CORE_UI_BRIEF_H_

#include "core/include/monte_carlo.h"
#include "core/include/criteria_metropolis.h"
#include "core/include/trial_translate.h"

namespace feasst {

inline void set_log(const std::string file_name,
    const int steps,
    MonteCarlo * mc) {
  auto log = std::make_shared<Log>();
  log->set_steps_per_write(steps);
  log->set_file_name(file_name);
  mc->add(log);
}

inline void set_energy_check(const double tolerance,
    const int steps,
    MonteCarlo * mc) {
  auto checker = std::make_shared<EnergyCheck>();
  checker->set_steps_per_update(steps);
  checker->set_tolerance(tolerance);
  mc->add(checker);
}

inline void set_trial_tune(const int steps, MonteCarlo * mc) {
  auto tuner = std::make_shared<Tuner>();
  tuner->set_steps_per_update(steps);
  mc->add(tuner);
}

inline void set_metropolis_criteria(
    const double beta,
    const double activity,
    MonteCarlo * mc) {
  auto criteria = std::make_shared<CriteriaMetropolis>();
  criteria->set_beta(beta);
  criteria->add_activity(activity);
  mc->set(criteria);
}

inline void set_trial_translate(const double weight,
    MonteCarlo * mc) {
  auto trial = std::make_shared<TrialTranslate>();
  trial->set_weight(weight);
  mc->add(trial);
}

}  // namespace feasst

#endif  // FEASST_CORE_UI_BRIEF_H_
