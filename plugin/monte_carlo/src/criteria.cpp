#include <cmath>
#include <sstream>
#include <iostream>
#include "utils/include/utils.h"  // is_equal
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/constraint.h"

namespace feasst {

Criteria::Criteria() {
  set_expanded_state();
  data_.get_dble_1D()->resize(1);
}

Criteria::Criteria(std::shared_ptr<Constraint> constraint) : Criteria() {
  add(constraint);
}

std::string Criteria::status_header() const {
  std::stringstream ss;
  if (num_states() > 1) {
    ss << ",state";
  }
  ss << ",energy";
  return ss.str();
}

std::string Criteria::status() const {
  std::stringstream ss;
  if (num_states() > 1) {
    ss << "," << state();
  }
  ss << "," << current_energy();
  return ss.str();
}

std::string Criteria::write() const {
  std::stringstream ss;
  return ss.str();
}

void Criteria::set_expanded_state(const int state, const int num) {
  expanded_state_ = state;
  num_expanded_states_ = num;
}

void Criteria::revert_(const bool accepted, const double ln_prob) {
  if (accepted) {
    *current_energy_() = previous_energy_;
  }
}


std::map<std::string, std::shared_ptr<Criteria> >& Criteria::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Criteria> >* ans =
     new std::map<std::string, std::shared_ptr<Criteria> >();
  return *ans;
}

void Criteria::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Criteria> Criteria::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Criteria> Criteria::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Criteria> Criteria::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Criteria> Criteria::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

bool Criteria::is_equal(const Criteria& criteria,
    const double tolerance) const {
  if (std::abs(current_energy() - criteria.current_energy()) > tolerance) {
    INFO(MAX_PRECISION << "current energy not equal: " << current_energy()
      << " vs " << criteria.current_energy() << " tol " << tolerance);
    return false;
  }
  if (expanded_state_ != criteria.expanded_state_) {
    INFO("expanded_states not equal: " << expanded_state_
      << " vs " << criteria.expanded_state_);
    return false;
  }
// HWH this doesn't consider tolerance
//  std::stringstream ss1, ss2;
//  serialize(ss1);
//  criteria->serialize(ss2);
//  if (ss1.str() != ss2.str()) {
//    INFO(ss1.str());
//    INFO(ss2.str());
//    return false;
//  }
  return true;
}

bool Criteria::is_equal(const Criteria& criteria) const {
  return is_equal(criteria, NEAR_ZERO);
}

void Criteria::serialize_criteria_(std::ostream& ostr) const {
  feasst_serialize_version(692, ostr);
  feasst_serialize(previous_energy_, ostr);
  feasst_serialize(expanded_state_, ostr);
  feasst_serialize(num_expanded_states_, ostr);
  feasst_serialize(phase_, ostr);
  feasst_serialize_fstdr(constraints_, ostr);
  feasst_serialize_fstobj(data_, ostr);
}

Criteria::Criteria(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 692, "version mismatch: " << version);
  feasst_deserialize(&previous_energy_, istr);
  feasst_deserialize(&expanded_state_, istr);
  feasst_deserialize(&num_expanded_states_, istr);
  feasst_deserialize(&phase_, istr);
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
  feasst_deserialize_fstobj(&data_, istr);
}

void Criteria::set_current_energy(const double energy) {
  previous_energy_ = current_energy();
  *current_energy_() = energy;
  DEBUG("setting current energy: " << current_energy());
  DEBUG("previous " << previous_energy_);
}

/// Return whether constraints are statisfied.
bool Criteria::is_allowed(const System& system, const Acceptance& acceptance) {
  for (const std::shared_ptr<Constraint> con : constraints_) {
    if (!con->is_allowed(system, *this, acceptance)) {
  //for (int con = 0; con < static_cast<int>(constraints_.size()); ++con) {
    //if (!constraints_[con]->is_allowed(system, this, acceptance)) {
      return false;
    }
  }
  return true;
}

void Criteria::set_num_iterations(const int iteration) {
  FATAL("not implemented");
}

}  // namespace feasst
