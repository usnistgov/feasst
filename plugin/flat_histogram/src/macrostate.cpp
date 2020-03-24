#include <cmath>
#include "utils/include/serialize.h"
#include "flat_histogram/include/macrostate.h"
#include "math/include/utils_math.h"

namespace feasst {

Macrostate::Macrostate(const Histogram& histogram, const argtype& args) {
  set(histogram);

  // soft limits
  args_.init(args);
  soft_min_ = 0;
  soft_max_ = histogram_.size() - 1;
  if (args_.key("soft_max").used()) {
    soft_max_ = args_.integer();
    if (args_.key("soft_min").used()) {
      soft_min_ = args_.integer();
    }
  }
  DEBUG("soft min " << soft_min_);
  DEBUG("soft max " << soft_max_);
  DEBUG("edges " << feasst_str(histogram_.edges()));
}

bool Macrostate::is_allowed(const System* system,
                            const Criteria* criteria,
                            const Acceptance& acceptance) {
  const double val = value(system, criteria);
  if (val > histogram_.max() || val < histogram_.min()) {
    return false;
  }
  const int ibin = histogram_.bin(val) + acceptance.macrostate_shift();
  DEBUG("ibin " << ibin << " max " << soft_max() << " min " << soft_min());
  if (ibin > soft_max() or ibin < soft_min()) {
    return false;
  }
  for (int con = 0; con < static_cast<int>(constraints_.size()); ++con) {
    if (!constraints_[con]->is_allowed(system, criteria, acceptance)) {
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
  feasst_serialize_fstdr(constraints_, ostr);
}

Macrostate::Macrostate(std::istream& istr) {
  ASSERT(feasst_deserialize_version(istr) == 520, "version check");
  feasst_deserialize_fstobj(&histogram_, istr);
  feasst_deserialize(&soft_max_, istr);
  feasst_deserialize(&soft_min_, istr);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize_fstdr(constraints_, istr);
  { int dim1;
    istr >> dim1;
    constraints_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstdr((*vector)[index], istr);
      int existing;
      istr >> existing;
      if (existing != 0) {
        constraints_[index] = constraints_[index]->deserialize(istr);
      }
    }
  }
}

std::vector<double> segment(
    const double min,
    const double max,
    const int num,
    const double exp) {
  ASSERT(num > 0, "num(" << num << ") must be > 0.");
  std::vector<double> segment(num + 1);
  segment[0] = min;
  segment[num] = max;
  const long double exp_diff =
    (std::pow(max, exp) - std::pow(min, exp))/static_cast<double>(num);
  for (int index = 1; index < num; ++index) {
    segment[index] = std::pow(std::pow(segment[index - 1], exp) + exp_diff, 1./exp);
  }
  return segment;
}

/// Segment a range into windows by exponential scaling.
std::vector<std::vector<int> > window(
    const int min,
    const int max,
    const int num,
    const double exp,
    const int extra_overlap) {
  std::vector<double> boundaries = segment(min, max, num, exp);
  std::vector<std::vector<int> > windows(num, std::vector<int>(2, 0.));
  for (int index = 0; index < num; ++index) {
    if (index == 0) {
      windows[index][0] = min;
    } else {
      windows[index][0] = round(boundaries[index] - extra_overlap);
    }
    windows[index][1] = round(boundaries[index + 1]);
  }
  return windows;
}

}  // namespace feasst
