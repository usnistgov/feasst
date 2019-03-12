
#include "core/include/debug.h"
#include "core/include/analyze.h"

namespace feasst {

void Analyze::trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  if (is_time(steps_per_update(), &steps_since_update_)) {
    update(criteria, system, trial_factory);
  }
  if (is_time(steps_per_write(), &steps_since_write_)) {
    printer(write(criteria, system, trial_factory));
  }
}

void Analyze::update(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ERROR("not implemented");
}

std::string Analyze::write(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ERROR("not implemented");
  return std::string("");
}

void AnalyzeFactory::initialize(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  for (const std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->initialize(criteria, system, trial_factory);
  }
}

void AnalyzeFactory::trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  for (const std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->trial(criteria, system, trial_factory);
  }
}

}  // namespace feasst
