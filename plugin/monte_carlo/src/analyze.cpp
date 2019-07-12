
#include "utils/include/debug.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Analyze> >& Analyze::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Analyze> >* ans =
     new std::map<std::string, std::shared_ptr<Analyze> >();
  return *ans;
}

void Analyze::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Analyze> Analyze::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Analyze> Analyze::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(),
    istr,
    // rewind reading of class name
    true);
}

void Analyze::trial(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  if (is_time(steps_per_update(), &steps_since_update_)) {
    update(criteria, system, trial_factory);
  }
  if (is_time(steps_per_write(), &steps_since_write_)) {
    printer(write(criteria, system, trial_factory));
  }
}

void Analyze::update(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ERROR("not implemented");
}

std::string Analyze::write(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ERROR("not implemented");
  return std::string("");
}

}  // namespace feasst
