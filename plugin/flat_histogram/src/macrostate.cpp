
#include "flat_histogram/include/macrostate.h"
#include "math/include/utils_math.h"

namespace feasst {

bool Macrostate::is_allowed(const System* system,
                            const Criteria* criteria,
                            const int shift) {
  const double val = value(system, criteria);
  if (val > histogram_.max() || val < histogram_.min()) {
    return false;
  }
  const int ibin = histogram_.bin(val) + shift;
  DEBUG("ibin " << ibin << " max " << soft_max() << " min " << soft_min());
  if (ibin > soft_max() or ibin < soft_min()) {
    return false;
  }
  for (int con = 0; con < static_cast<int>(constraints_.size()); ++con) {
    if (!constraints_[con]->is_allowed(system, criteria)) {
      return false;
    }
  }
  return true;
}

void Macrostate::swap_soft_bounds(Macrostate * macrostate) {
  swap(&soft_max_, &(macrostate->soft_max_));
  swap(&soft_min_, &(macrostate->soft_min_));
}

std::map<std::string, std::shared_ptr<Macrostate> >& Macrostate::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Macrostate> >* ans =
     new std::map<std::string, std::shared_ptr<Macrostate> >();
  return *ans;
}

void Macrostate::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Macrostate> Macrostate::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Macrostate> Macrostate::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void Macrostate::serialize_macrostate_(std::ostream& ostr) const {
  feasst_serialize_version(520, ostr);
  feasst_serialize_fstobj(histogram_, ostr);
  feasst_serialize(soft_max_, ostr);
  feasst_serialize(soft_min_, ostr);
  ASSERT(constraints_.size() == 0, "constraint serialization not implemented");
}

Macrostate::Macrostate(std::istream& istr) {
  ASSERT(feasst_deserialize_version(istr) == 520, "version check");
  feasst_deserialize_fstobj(&histogram_, istr);
  feasst_deserialize(&soft_max_, istr);
  feasst_deserialize(&soft_min_, istr);
}

}  // namespace feasst
