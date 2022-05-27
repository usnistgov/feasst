#include "utils/include/serialize.h"
#include "monte_carlo/include/constraint.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Constraint> >&
  Constraint::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Constraint> >* ans =
     new std::map<std::string, std::shared_ptr<Constraint> >();
  return *ans;
}

void Constraint::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Constraint> Constraint::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Constraint> Constraint::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Constraint> Constraint::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

std::shared_ptr<Constraint> Constraint::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void Constraint::serialize_constraint_(std::ostream& ostr) const {
  feasst_serialize_version(4176, ostr);
}

Constraint::Constraint(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(4176 == version, "mismatch version: " << version);
}

}  // namespace feasst

