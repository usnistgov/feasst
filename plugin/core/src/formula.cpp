#include "core/include/formula.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Formula> >& Formula::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Formula> >* ans =
     new std::map<std::string, std::shared_ptr<Formula> >();
  return *ans;
}

void Formula::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Formula> Formula::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Formula> Formula::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void Formula::serialize_formula_(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
  feasst_serialize(x0_, ostr);
}

Formula::Formula(std::istream& istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&x0_, istr);
}

}  // namespace feasst
