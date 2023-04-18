#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Modify> >& Modify::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Modify> >* ans =
     new std::map<std::string, std::shared_ptr<Modify> >();
  return *ans;
}

void Modify::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Modify> Modify::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Modify> Modify::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Modify> Modify::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(),
    istr,
    // rewind reading of class name
    true);
}

std::shared_ptr<Modify> Modify::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

void Modify::check_update_(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  DEBUG("check update " << trials_per_update() << " " << trials_since_update());
  if (is_time(trials_per_update(), &trials_since_update_)) {
    update(criteria, system, trial_factory);
  }
}

void Modify::trial(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  if ((stop_after_phase() == -1 || criteria->phase() <= stop_after_phase()) &&
      (stop_after_iteration() == -1 || criteria->num_iterations() <= stop_after_iteration())) {
    if ((criteria->phase() > start_after_phase()) &&
        (criteria->num_iterations() > start_after_iteration())) {
      check_update_(criteria, system, trial_factory);
      if (is_time(trials_per_write(), &trials_since_write_)) {
        write_to_file(criteria, system, trial_factory);
      }
    }
  }
}

void Modify::write_to_file(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  if (trials_per_write() != -1) {
    printer(write(criteria, system, trial_factory), file_name(*criteria));
  }
}

void Modify::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  FATAL("not implemented");
}

std::string Modify::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  DEBUG(trials_per_write());
  FATAL(class_name() << " not implemented");
  return std::string("");
}

const std::vector<std::shared_ptr<Modify> >& Modify::modifiers() const {
  FATAL("not implemented");
}

const Modify& Modify::modify(const int index) const {
  FATAL("not implemented");
}

Modify * Modify::get_modify(const int index) {
  FATAL("not implemented");
}

ModifyUpdateOnly::ModifyUpdateOnly(argtype * args) : Modify(args) {
  // disable write
  Modify::set_trials_per_write(-1);

  // parse
  if (used("trials_per", *args)) {
    WARN("ModifyUpdateOnly::trials_per is deprecated. Use trials_per_update.");
    set_trials_per(integer("trials_per", args));
  }
}

void ModifyUpdateOnly::set_trials_per_write(const int trials) {
  ERROR("This modify is update only.");
}

}  // namespace feasst
