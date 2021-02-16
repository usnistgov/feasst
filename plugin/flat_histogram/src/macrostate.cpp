#include <cmath>
#include "utils/include/serialize.h"
#include "flat_histogram/include/macrostate.h"
#include "math/include/utils_math.h"

namespace feasst {

Macrostate::Macrostate(const Histogram& histogram, argtype args)
  : Macrostate(histogram, &args) {
  check_all_used(args);
}
Macrostate::Macrostate(const Histogram& histogram, argtype * args) {
  set(histogram);

  // soft limits
  soft_min_ = 0;
  soft_max_ = histogram_.size() - 1;
  if (used("soft_max", *args)) {
    soft_max_ = integer("soft_max", args);
    if (used("soft_min", *args)) {
      soft_min_ = integer("soft_min", args);
    }
  }
  DEBUG("soft min " << soft_min_);
  DEBUG("soft max " << soft_max_);
  DEBUG("edges " << feasst_str(histogram_.edges()));
}

bool Macrostate::is_allowed(const System& system,
                            const Criteria& criteria,
                            const Acceptance& acceptance) const {
  const double val = value(system, criteria, acceptance);
  if (val > histogram_.max() || val < histogram_.min()) {
    return false;
  }
  const int ibin = histogram_.bin(val);// + acceptance.macrostate_shift();
  DEBUG("ibin " << ibin << " max " << soft_max() << " min " << soft_min());
  if (ibin > soft_max() or ibin < soft_min()) {
    return false;
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
}

Macrostate::Macrostate(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 520, "version: " << version);
  feasst_deserialize_fstobj(&histogram_, istr);
  feasst_deserialize(&soft_max_, istr);
  feasst_deserialize(&soft_min_, istr);
}

}  // namespace feasst
