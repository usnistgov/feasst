#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Analyze> >& Analyze::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Analyze> >* ans =
     new std::map<std::string, std::shared_ptr<Analyze> >();
  return *ans;
}

void Analyze::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Analyze> Analyze::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Analyze> Analyze::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Analyze> Analyze::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(),
    istr,
    // rewind reading of class name
    true);
}

std::shared_ptr<Analyze> Analyze::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

void Analyze::check_update_(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  if (is_time(steps_per_update(), &steps_since_update_)) {
    update(criteria, system, trial_factory);
  }
}

void Analyze::trial(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  if (stop_after_phase() == -1 ||
      criteria.phase() <= stop_after_phase()) {
    if (criteria.phase() > start_after_phase()) {
      check_update_(criteria, system, trial_factory);
      if (is_time(steps_per_write(), &steps_since_write_)) {
        printer(write(criteria, system, trial_factory), file_name(criteria));
      }
    }
  }
}

void Analyze::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  FATAL("not implemented");
}

std::string Analyze::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  FATAL("not implemented");
  return std::string("");
}

const std::vector<std::shared_ptr<Analyze> >& Analyze::analyzers() const {
  FATAL("not implemented");
}

const Analyze& Analyze::analyze(const int index) const {
  FATAL("not implemented");
}

AnalyzeWriteOnly::AnalyzeWriteOnly(argtype * args) : Analyze(args) {
  // disable update
  Stepper::set_steps_per_update(-1);

  // parse
  if (used("steps_per", *args)) {
    set_steps_per(integer("steps_per", args));
  }
}

void AnalyzeWriteOnly::set_steps_per_update(const int steps) {
  ERROR("This analyze is write only.");
}

AnalyzeUpdateOnly::AnalyzeUpdateOnly(argtype * args) : Analyze(args) {
  ASSERT(steps_per_write() == 1, "AnalyzeUpdateOnly doesn't the argument " <<
    "steps_per_write");
  ASSERT(file_name().empty(), "AnalyzeUpdateOnly doesn't the argument " <<
    "file_name");

  // disable write
  Analyze::set_steps_per_write(-1);

  // parse
  if (used("steps_per", *args)) {
    set_steps_per(integer("steps_per", args));
  }
}

void AnalyzeUpdateOnly::set_steps_per_write(const int steps) {
  ERROR("This analyze is update only.");
}

}  // namespace feasst
