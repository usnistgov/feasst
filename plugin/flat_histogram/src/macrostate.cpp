#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize_extra.h"
#include "math/include/utils_math.h"
#include "math/include/histogram.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/acceptance.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/macrostate_energy.h"

namespace feasst {

Macrostate::Macrostate(const Histogram& histogram, argtype args)
  : Macrostate(histogram, &args) {
  feasst_check_all_used(args);
}
Macrostate::Macrostate(const Histogram& histogram, argtype * args) {
  set(histogram);
  // soft limits
  soft_min_ = 0;
  soft_max_ = histogram_->size() - 1;
  if (used("soft_macro_max", *args)) {
    soft_max_ = histogram_->bin(dble("soft_macro_max", args));
    if (used("soft_macro_min", *args)) {
      soft_min_ = histogram_->bin(dble("soft_macro_min", args));
    }
  }
  DEBUG("soft min " << soft_min_);
  DEBUG("soft max " << soft_max_);
  DEBUG("edges " << feasst_str(histogram_->edges()));
}

Macrostate::Macrostate(argtype args) :
    Macrostate(Histogram(&args), &args) {
  feasst_check_all_used(args);
}

bool Macrostate::is_allowed(const System& system,
                            const Criteria& criteria,
                            const Acceptance& acceptance) const {
  const double val = value(system, criteria, acceptance);
  if (val > histogram_->max() || val < histogram_->min()) {
    return false;
  }
  const int ibin = histogram_->bin(val);// + acceptance.macrostate_shift();
  DEBUG("ibin " << ibin << " max " << soft_max() << " min " << soft_min());
  if (ibin > soft_max() or ibin < soft_min()) {
    return false;
  }
  return true;
}

//void Macrostate::swap_soft_bounds(Macrostate * macrostate) {
//  feasst_swap(&soft_max_, &(macrostate->soft_max_));
//  feasst_swap(&soft_min_, &(macrostate->soft_min_));
//}

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
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void Macrostate::serialize_macrostate_(std::ostream& ostr) const {
  feasst_serialize_version(520, ostr);
  feasst_serialize(histogram_, ostr);
  feasst_serialize(soft_max_, ostr);
  feasst_serialize(soft_min_, ostr);
}

Macrostate::Macrostate(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 520, "version: " << version);
  feasst_deserialize(histogram_, istr);
  feasst_deserialize(&soft_max_, istr);
  feasst_deserialize(&soft_min_, istr);
}

std::shared_ptr<Macrostate> Macrostate::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

std::shared_ptr<Macrostate> Macrostate::create(argtype * args) const {
  FATAL("not implemented");
}

int Macrostate::set_soft_max(const int index, const System& sys, const Criteria& criteria) {
  Acceptance acc;
  if (bin(sys, criteria, acc) <= index) {
    soft_max_ = index;
    return 1;
  } else {
    return 0;
  }
}

int Macrostate::set_soft_min(const int index, const System& sys, const Criteria& criteria) {
  Acceptance acc;
  if (bin(sys, criteria, acc) >= index) {
    soft_min_ = index;
    return 1;
  } else {
    return 0;
  }
}

double Macrostate::value(const int bin) const {
  return histogram_->center_of_bin(bin + soft_min_);
}

void Macrostate::set(const Histogram histogram) { histogram_ = std::make_unique<Histogram>(histogram); }
const Histogram& Macrostate::histogram() const { return *histogram_; }
int Macrostate::bin(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  return histogram_->bin(value(system, criteria, acceptance));
}

}  // namespace feasst
