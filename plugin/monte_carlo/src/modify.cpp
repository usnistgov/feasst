#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

Modify::~Modify() {}

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

void Modify::check_update_(MonteCarlo * mc) {
  DEBUG("check update " << trials_per_update() << " " << trials_since_update());
  if (is_time(trials_per_update(), &trials_since_update_)) {
    update(mc);
  }
}

void Modify::trial(MonteCarlo * mc) {
  const Criteria& criteria = mc->criteria();
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

void Modify::write_to_file(MonteCarlo * mc) {
  if (trials_per_write() != -1) {
    printer(write(mc), output_file(mc->criteria()));
  }
}

void Modify::update(MonteCarlo * mc) {
  FATAL("not implemented");
}

std::string Modify::write(MonteCarlo * mc) {
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
  ASSERT(output_file().empty(),
    "ModifyUpdateOnly does not use the argument output_file.");
}

void ModifyUpdateOnly::set_trials_per_write(const int trials) {
  ERROR("This modify is update only.");
}

}  // namespace feasst
