#include "utils/include/serialize_extra.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/tunable.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/perturb.h"

namespace feasst {

Perturb::Perturb(argtype args) : Perturb(&args) { feasst_check_all_used(args); }
Perturb::Perturb(argtype * args) {
  tunable_ = std::make_shared<Tunable>(args);
}

void Perturb::before_select() {
  revert_possible_ = false;
  finalize_possible_  = false;
}

void Perturb::set_revert_possible(const bool revert_possible,
    TrialSelect * revert_select) {
  revert_possible_ = revert_possible;
  revert_select_ = revert_select;
}

void Perturb::revert(System * system) {
  if (revert_possible_) {
    FATAL("not implemented");
  }
}

void Perturb::set_finalize_possible(const bool finalize_possible,
  /// If possible, store the selection.
  TrialSelect * finalize_select) {
  finalize_possible_ = finalize_possible;
  finalize_select_ = finalize_select;
  DEBUG("finalize possible: " << finalize_possible_);
}

void Perturb::finalize(System * system) {
  //system->finalize();
  if (finalize_possible_) {
    FATAL("not implemented");
  }
}

std::string Perturb::status_header() const {
  std::stringstream ss;
  if (tunable().is_enabled()) {
    ss << ",tunable";
  }
  return ss.str();
}

std::string Perturb::status() const {
  std::stringstream ss;
  if (tunable().is_enabled()) {
    ss << "," << tunable_->value();
  }
  return ss.str();
}

std::map<std::string, std::shared_ptr<Perturb> >& Perturb::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Perturb> >* ans =
     new std::map<std::string, std::shared_ptr<Perturb> >();
  return *ans;
}

void Perturb::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Perturb> Perturb::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Perturb> Perturb::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void Perturb::serialize_perturb_(std::ostream& ostr) const {
  feasst_serialize_version(902, ostr);
  feasst_serialize(tunable_, ostr);
}

Perturb::Perturb(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(902 == version, "mismatch version: " << version);
//  feasst_deserialize(tunable_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      tunable_ = std::make_shared<Tunable>(istr);
    }
  }
}

void Perturb::perturb(
  System * system,
  TrialSelect * select,
  Random * random,
  const bool is_position_held,
  Acceptance * acceptance) {
  FATAL("not implemented");
}

const Tunable& Perturb::tunable() const { return *tunable_; }
void Perturb::set_tunable(const double value) { tunable_->set_value(value); }
void Perturb::tune(const double actual) { tunable_->tune(actual); }
void Perturb::disable_tunable_() { tunable_->disable(); }

void Perturb::set_tunable_min_and_max(const double min, const double max) {
  tunable_->set_min_and_max(min, max);
}

}  // namespace feasst
