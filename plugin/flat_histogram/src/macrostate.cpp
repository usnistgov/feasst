
#include "flat_histogram/include/macrostate.h"

namespace feasst {

bool Macrostate::is_in_range(const System* system, const Criteria* criteria) {
  const double val = value(system, criteria);
  if (val <= histogram_.max() &&
      val >= histogram_.min()) {
    return true;
  }
  return false;
}

std::map<std::string, std::shared_ptr<Macrostate> >& Macrostate::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Macrostate> >* ans =
     new std::map<std::string, std::shared_ptr<Macrostate> >();
  return *ans;
}

void Macrostate::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Macrostate> Macrostate::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Macrostate> Macrostate::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void Macrostate::serialize_macrostate_(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
  feasst_serialize_fstobj(histogram_, ostr);
}

Macrostate::Macrostate(std::istream& istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize_fstobj(&histogram_, istr);
}

}  // namespace feasst
