
#include "utils/include/debug.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Modify> >& Modify::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Modify> >* ans =
     new std::map<std::string, std::shared_ptr<Modify> >();
  return *ans;
}

void Modify::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Modify> Modify::create(std::istream& istr) const {
  ERROR("not implemented");
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
  if (is_time(steps_per_update(), &steps_since_update_)) {
    update(criteria, system, trial_factory);
  }
  if (is_time(steps_per_write(), &steps_since_write_)) {
    printer(write(criteria, system, trial_factory));
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
}  // namespace feasst
