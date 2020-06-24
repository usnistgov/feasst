#include "math/include/formula.h"
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"
#include "utils/include/serialize.h"

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

std::shared_ptr<Formula> Formula::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

Formula::Formula(const argtype& args) {
  args_.init(args);
  x0_ = args_.key("x0").dflt("0").dble();
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
