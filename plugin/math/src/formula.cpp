#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/arguments.h"
#include "math/include/formula.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Formula> >& Formula::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Formula> >* ans =
     new std::map<std::string, std::shared_ptr<Formula> >();
  return *ans;
}

void Formula::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Formula> Formula::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Formula> Formula::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Formula> Formula::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Formula> Formula::factory(const std::string name,
                                          argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

Formula::Formula(argtype * args) {
  x0_ = dble("x0", args, 0.);
}

void Formula::serialize_formula_(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(5694, ostr);
  feasst_serialize(x0_, ostr);
}

Formula::Formula(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5694, "version mismatch: " << version);
  feasst_deserialize(&x0_, istr);
}

double Formula::evaluate(const double x) const { FATAL("not implemented"); }

double Formula::evaluate(const double x, const double y) const {
  FATAL("not implemented"); }

double Formula::derivative(const double x) const { FATAL("not implemented"); }

}  // namespace feasst
