
#include "flat_histogram/include/macrostate.h"
#include "math/include/utils_math.h"

namespace feasst {

bool Macrostate::is_in_range(const System* system, const Criteria* criteria) {
  const double val = value(system, criteria);
  if (val <= soft_max() and val >= soft_min()) {
    return true;
  } else {
    return false;
  }
}

const int Macrostate::soft_max() const {
  if (is_soft_bound_) {
    return soft_max_;
  } else {
    return histogram_.max();
  }
}

const int Macrostate::soft_min() const {
  if (is_soft_bound_) {
    return soft_min_;
  } else {
    return histogram_.min();
  }
}

void Macrostate::swap_soft_bounds(Macrostate * macrostate) {
  ASSERT(is_soft_bound_ and macrostate->is_soft_bound_,
    "macrostates are not both soft bound");
  swap(&soft_max_, &(macrostate->soft_max_));
  swap(&soft_min_, &(macrostate->soft_min_));
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
  feasst_serialize_version(520, ostr);
  feasst_serialize_fstobj(histogram_, ostr);
  feasst_serialize(is_soft_bound_, ostr);
  feasst_serialize(soft_max_, ostr);
  feasst_serialize(soft_min_, ostr);
}

Macrostate::Macrostate(std::istream& istr) {
  ASSERT(feasst_deserialize_version(istr) == 520, "version check");
  feasst_deserialize_fstobj(&histogram_, istr);
  feasst_deserialize(&is_soft_bound_, istr);
  feasst_deserialize(&soft_max_, istr);
  feasst_deserialize(&soft_min_, istr);
}

}  // namespace feasst
