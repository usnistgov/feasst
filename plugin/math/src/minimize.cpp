#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/minimize.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Minimize> >& Minimize::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Minimize> >* ans =
     new std::map<std::string, std::shared_ptr<Minimize> >();
  return *ans;
}

void Minimize::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Minimize> Minimize::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Minimize> Minimize::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

Minimize::Minimize(argtype args) {
  tolerance_ = dble("tolerance", &args);
  set_lower(dble("lower", &args));
  set_upper(dble("upper", &args));
  FEASST_CHECK_ALL_USED(args);
}

void Minimize::serialize_solver_(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(6236, ostr);
  feasst_serialize(tolerance_, ostr);
  feasst_serialize(lower_, ostr);
  feasst_serialize(upper_, ostr);
}

Minimize::Minimize(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6236, "version mismatch: " << version);
  feasst_deserialize(&tolerance_, istr);
  feasst_deserialize(&lower_, istr);
  feasst_deserialize(&upper_, istr);
}

double Minimize::minimum(Formula * formula) {
  double a, b;
  bracket(&a, &b, formula);
  DEBUG("a: " << a << " b: " << b);
  return 0.5*(a + b);
}

}  // namespace feasst
