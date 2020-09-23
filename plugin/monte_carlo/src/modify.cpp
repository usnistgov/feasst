#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

Modify::Modify(const argtype &args) : Stepper(args) {
  ASSERT(!is_multistate() || !is_multistate_aggregate(),
    "multistate_aggregate not implemented");
}

std::map<std::string, std::shared_ptr<Modify> >& Modify::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Modify> >* ans =
     new std::map<std::string, std::shared_ptr<Modify> >();
  return *ans;
}

void Modify::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Modify> Modify::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Modify> Modify::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(),
    istr,
    // rewind reading of class name
    true);
}

void Modify::trial(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  if (stop_after_phase() == -1 ||
      criteria->phase() <= stop_after_phase()) {
    if (criteria->phase() > start_after_phase()) {
      if (is_time(steps_per_update(), &steps_since_update_)) {
        update(criteria, system, trial_factory);
      }
      if (is_time(steps_per_write(), &steps_since_write_)) {
        printer(write(criteria, system, trial_factory),
                file_name(static_cast<const Criteria&>(*criteria)));
      }
    }
  }
}

ModifyUpdateOnly::ModifyUpdateOnly(const argtype &args) : Modify(args) {
  // disable write
  Modify::set_steps_per_write(-1);

  // parse
  if (!args_.key("steps_per").empty()) {
    set_steps_per(args_.integer());
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
  FATAL("not implemented");
  return std::string("");
}

const std::vector<std::shared_ptr<Modify> >& Modify::modifiers() const {
  FATAL("not implemented");
}

const Modify& Modify::modify(const int index) const {
  FATAL("not implemented");
}

void ModifyUpdateOnly::set_steps_per_write(const int steps) {
  ERROR("This modify is update only.");
}

}  // namespace feasst
