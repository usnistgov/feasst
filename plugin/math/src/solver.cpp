#include "math/include/solver.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Solver> >& Solver::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Solver> >* ans =
     new std::map<std::string, std::shared_ptr<Solver> >();
  return *ans;
}

void Solver::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Solver> Solver::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Solver> Solver::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

Solver::Solver(argtype args) {
  tolerance_ = dble("tolerance", &args);
  if (used("lower", args)) set_lower(dble("lower", &args));
  if (used("upper", args)) set_upper(dble("upper", &args));
  if (used("guess", args)) set_guess(dble("guess", &args));
}

void Solver::serialize_solver_(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(3756, ostr);
  feasst_serialize(tolerance_, ostr);
  feasst_serialize(lower_, ostr);
  feasst_serialize(is_lower_, ostr);
  feasst_serialize(upper_, ostr);
  feasst_serialize(is_upper_, ostr);
  feasst_serialize(guess_, ostr);
  feasst_serialize(is_guess_, ostr);
}

Solver::Solver(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3756, "version mismatch: " << version);
  feasst_deserialize(&tolerance_, istr);
  feasst_deserialize(&lower_, istr);
  feasst_deserialize(&is_lower_, istr);
  feasst_deserialize(&upper_, istr);
  feasst_deserialize(&is_upper_, istr);
  feasst_deserialize(&guess_, istr);
  feasst_deserialize(&is_guess_, istr);
}

double Solver::root(Formula * formula) {
  FATAL("not implemented");
}

}  // namespace feasst
