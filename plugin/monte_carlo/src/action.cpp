#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "monte_carlo/include/action.h"

namespace feasst {

Action::Action(argtype * args) {}
Action::Action(argtype args) : Action(&args) { check_all_used(args); }

std::map<std::string, std::shared_ptr<Action> >& Action::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Action> >* ans =
     new std::map<std::string, std::shared_ptr<Action> >();
  return *ans;
}

void Action::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Action> Action::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Action> Action::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Action> Action::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Action> Action::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

void Action::serialize_action_(std::ostream& ostr) const {
  feasst_serialize_version(6937, ostr);
}

Action::Action(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(6937 == version, "mismatch version: " << version);
}

void Action::perform(MonteCarlo * mc) { FATAL("not implemented"); }

}  // namespace feasst
