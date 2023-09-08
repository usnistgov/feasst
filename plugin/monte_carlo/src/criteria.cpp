#include <cmath>
#include <sstream>
#include <iostream>
#include "utils/include/utils.h"  // is_equal
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/constraint.h"
#include "monte_carlo/include/constrain_num_particles.h"

namespace feasst {

Criteria::Criteria(argtype * args) {
  set_expanded_state();
  data_.get_dble_2D()->resize(1);
  *current_energy_() = std::vector<double>();
  data_.get_dble_3D()->resize(1);
  *current_energy_profile_() = std::vector<std::vector<double> >();
  data_.get_int_1D()->resize(2);
  *num_iterations_() = 0;
  *num_attempt_since_last_iteration_() = 0;
  num_iterations_to_complete_ = integer("num_iterations_to_complete", args, 20);
  if (used("Constraint", *args)) {
    add(ConstrainNumParticles().factory(str("Constraint", args), args));
  }
}
Criteria::Criteria(argtype args) : Criteria(&args) {
  FEASST_CHECK_ALL_USED(args);
}

Criteria::Criteria(std::shared_ptr<Constraint> constraint, argtype args)
  : Criteria(args) {
  add(constraint);
}

double Criteria::current_energy(const int config) const {
  ASSERT(config < static_cast<int>(data_.dble_2D()[0].size()),
    "config:" << config << " >= size:" << data_.dble_2D()[0].size());
  return data_.dble_2D()[0][config];
}

const std::vector<double>& Criteria::current_energy_profile(const int config) const {
  ASSERT(config < static_cast<int>(data_.dble_3D()[0].size()),
    "config:" << config << " >= size:" << data_.dble_3D()[0].size());
//  if (config == static_cast<int>(const_current_energy_profile_()->size())) {
//    const_current_energy_profile_()->resize(config + 1);
//  }
  return data_.dble_3D()[0][config];
}

void Criteria::update_current_energy(const Acceptance& acceptance) {
  for (int iconf = 0; iconf < acceptance.num_configurations(); ++iconf) {
    if (acceptance.updated(iconf) == 1) {
      DEBUG("iconf " << iconf);
      set_current_energy(acceptance.energy_new(iconf), iconf);
      set_current_energy_profile(acceptance.energy_profile_new(iconf), iconf);
    }
  }
}

std::string Criteria::status_header(const System& system) const {
  std::stringstream ss;
  if (num_states() > 1) {
    ss << ",state";
  }
  std::string append = "";
  for (int iconf = 0; iconf < system.num_configurations(); ++iconf) {
    if (system.num_configurations() > 1) {
      append = "_config" + str(iconf);
    }
    ss << ",energy" << append;
    for (int i = 0; i < static_cast<int>(current_energy_profile().size()); ++i) {
      std::string name;
      name = system.potential(i).model().class_name();
      if (name == "ModelEmpty") {
        name = system.potential(i).visit_model().class_name();
      }
      ss << "," << name << append;
    }
  }
  return ss.str();
}

std::string Criteria::status(const bool max_precision) const {
  std::stringstream ss;
  if (max_precision) {
    ss << MAX_PRECISION;
  }
  if (num_states() > 1) {
    ss << "," << state();
  }
  const int num_configs = static_cast<int>(const_current_energy_().size());
  for (int iconf = 0; iconf < num_configs; ++iconf) {
    ss << "," << current_energy(iconf);
    for (const double potential : current_energy_profile(iconf)) {
      ss << "," << potential;
    }
  }
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

void Criteria::revert_(const bool accepted, const bool endpoint, const double ln_prob) {
  if (accepted) {
    (*current_energy_())[0] = previous_energy_;
    (*current_energy_profile_())[0] = previous_energy_profile_;
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
  feasst_serialize(previous_energy_profile_, ostr);
  feasst_serialize(expanded_state_, ostr);
  feasst_serialize(num_expanded_states_, ostr);
  feasst_serialize(num_iterations_to_complete_, ostr);
  feasst_serialize(phase_, ostr);
  feasst_serialize_fstdr(constraints_, ostr);
  feasst_serialize_fstobj(data_, ostr);
}

Criteria::Criteria(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 692, "version mismatch: " << version);
  feasst_deserialize(&previous_energy_, istr);
  feasst_deserialize(&previous_energy_profile_, istr);
  feasst_deserialize(&expanded_state_, istr);
  feasst_deserialize(&num_expanded_states_, istr);
  feasst_deserialize(&num_iterations_to_complete_, istr);
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

void Criteria::set_current_energy(const double energy, const int config) {
  if (config == static_cast<int>(current_energy_()->size())) {
    current_energy_()->resize(config + 1);
  }
  previous_energy_ = current_energy(config);
  (*current_energy_())[config] = energy;
  DEBUG("setting current energy: " << current_energy(config));
  DEBUG("previous " << previous_energy_);
}

void Criteria::set_current_energy_profile(const std::vector<double>& energy,
    const int config) {
  if (config == static_cast<int>(current_energy_profile_()->size())) {
    current_energy_profile_()->resize(config + 1);
  }
  previous_energy_profile_ = current_energy_profile(config);
  (*current_energy_profile_())[config] = energy;
}

/// Return whether constraints are statisfied.
bool Criteria::is_allowed(const System& system, const Acceptance& acceptance) {
  for (const std::shared_ptr<Constraint>& con : constraints_) {
    if (!con->is_allowed(system, *this, acceptance)) {
  //for (int con = 0; con < static_cast<int>(constraints_.size()); ++con) {
    //if (!constraints_[con]->is_allowed(system, this, acceptance)) {
      return false;
    }
  }
  return true;
}

void Criteria::check_num_iterations_(const int num_trials_per_iteration) {
  *num_attempt_since_last_iteration_() += 1;
  if (*num_attempt_since_last_iteration_() >= num_trials_per_iteration) {
    *num_iterations_() += 1;
    *num_attempt_since_last_iteration_() = 0;
  }
}

int Criteria::set_soft_max(const int index, const System& sys) {
  FATAL("not implemented");
}

int Criteria::set_soft_min(const int index, const System& sys) {
  FATAL("not implemented");
}

void Criteria::set_cm(const bool inc_max, const int macro, const Criteria& crit) {
  FATAL("not implemented");
}

void Criteria::adjust_bounds(const bool left_most, const bool right_most,
  const bool left_complete, const bool right_complete,
  const bool all_min_size,
  const int min_size, const System& system, const System * upper_sys,
  Criteria * criteria, bool * adjusted_up, std::vector<int> * states) {
  FATAL("not implemented");
}

const Macrostate& Criteria::macrostate() const { FATAL("not implemented"); }
int Criteria::soft_min() const { FATAL("not implemented"); }
int Criteria::soft_max() const { FATAL("not implemented"); }

const Bias& Criteria::bias() const { FATAL("not implemented"); }

const FlatHistogram& Criteria::flat_histogram() const { FATAL("not implemented"); }

int Criteria::num_iterations(const int state) const {
  return const_num_iterations_();
}

}  // namespace feasst
