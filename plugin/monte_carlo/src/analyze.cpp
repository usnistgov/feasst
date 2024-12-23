#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
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

void Analyze::check_update_(const MonteCarlo& mc) {
  if (is_time(trials_per_update(), &trials_since_update_)) {
    update(mc);
  }
}

void Analyze::trial(const MonteCarlo& mc) {
  const Criteria& criteria = mc.criteria();
  if ((stop_after_phase() == -1 || criteria.phase() <= stop_after_phase()) &&
      (stop_after_cycle() == -1 || criteria.num_cycles() <= stop_after_cycle())) {
    if ((criteria.phase() > start_after_phase()) &&
        (criteria.num_cycles() > start_after_cycle())) {
      check_update_(mc);
      if (is_time(trials_per_write(), &trials_since_write_)) {
        write_to_file(mc);
      }
    }
  }
}

void Analyze::write_to_file(const MonteCarlo& mc) {
  if (trials_per_write() != -1) {
    printer(write(mc), output_file(mc.criteria()));
  }
}

void Analyze::update(const MonteCarlo& mc) {
  FATAL("not implemented");
}

std::string Analyze::write(const MonteCarlo& mc) {
  FATAL("not implemented");
  return std::string("");
}

const std::vector<std::shared_ptr<Analyze> >& Analyze::analyzers() const {
  FATAL("not implemented");
}

const Analyze& Analyze::analyze(const int index) const {
  FATAL("not implemented");
}

Analyze * Analyze::get_analyze(const int index) {
  FATAL("not implemented");
}

AnalyzeWriteOnly::AnalyzeWriteOnly(argtype * args) : Analyze(args) {
  // disable update
  Stepper::set_trials_per_update(-1);

  // parse
  if (used("trials_per", *args)) {
    WARN("AnalyzeWriteOnly::trials_per is deprecated. Use trials_per_write.");
    set_trials_per(integer("trials_per", args));
  }
}

void AnalyzeWriteOnly::set_trials_per_update(const int trials) {
  ERROR("This analyze is write only.");
}

AnalyzeUpdateOnly::AnalyzeUpdateOnly(argtype * args) : Analyze(args) {
  ASSERT(trials_per_write() == 1,
    "AnalyzeUpdateOnly does not use the argument trials_per_write.");
  ASSERT(output_file().empty(),
    "AnalyzeUpdateOnly does not use the argument output_file.");

  // disable write
  Analyze::set_trials_per_write(-1);

  // parse
  if (used("trials_per", *args)) {
    set_trials_per(integer("trials_per", args));
  }
}

void AnalyzeUpdateOnly::set_trials_per_write(const int trials) {
  ERROR("This analyze is update only.");
}

}  // namespace feasst
