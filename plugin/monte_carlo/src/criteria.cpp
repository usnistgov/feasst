#include "utils/include/arguments.h"
#include "utils/include/utils.h"  // is_equal
#include "utils/include/serialize_extra.h"
#include "utils/include/io.h"
#include "math/include/constants.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model.h"
#include "system/include/potential.h"
#include "system/include/model.h"
#include "system/include/visit_model_inner.h"
#include "system/include/system.h"
#include "system/include/bond_visitor.h"
#include "monte_carlo/include/acceptance.h"
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
  *num_cycles_() = 0;
  *num_attempt_since_last_cycle_() = 0;
  if (used("num_iterations_to_complete", *args)) {
    WARN("num_iterations_to_complete is deprecated. " <<
         "use cycles_to_complete.");
    cycles_to_complete_ = integer("num_iterations_to_complete", args);
  } else {
    cycles_to_complete_ = integer("cycles_to_complete", args, 20);
  }
  if (used("Constraint", *args)) {
    add(ConstrainNumParticles().factory(str("Constraint", args), args));
  }
}
Criteria::Criteria(argtype args) : Criteria(&args) {
  feasst_check_all_used(args);
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

std::string Criteria::status_header(const System& system,
    const bool include_bonds) const {
  std::stringstream ss;
  if (num_states() > 1) {
    ss << ",state";
  }
  std::string append = "";
  for (int iconf = 0; iconf < system.num_configurations(); ++iconf) {
    const std::string cname = system.configuration(iconf).name();
    if (system.num_configurations() > 1) {
      append = "_" + cname;
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
    if (include_bonds) {
      std::string append = "";
      if (system.num_configurations() > 1) {
        append = "_" + cname;
      }
      ss << ",BondTwoBody" << append <<
            ",BondThreeBody" << append <<
            ",BondFourBody" << append;
    }
  }
  return ss.str();
}

std::string Criteria::status(const System& system, const bool max_precision,
    const bool include_bonds, BondVisitor * visitor) const {
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
    if (include_bonds) {
      visitor->compute_all(system.configuration(iconf));
      if (max_precision) {
        ss << "," << MAX_PRECISION << visitor->energy_two_body()
           << "," << MAX_PRECISION << visitor->energy_three_body()
           << "," << MAX_PRECISION << visitor->energy_four_body();
      } else {
        ss << "," << visitor->energy_two_body()
           << "," << visitor->energy_three_body()
           << "," << visitor->energy_four_body();
      }
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

void Criteria::revert_(const bool accepted, const bool endpoint, const double ln_prob, const std::vector<int>& updated) {
  if (accepted) {
    for (int conf = 0; conf < static_cast<int>(updated.size()); ++conf) {
      if (updated[conf] == 1) {
        (*current_energy_())[conf] = previous_energy_[conf];
        (*current_energy_profile_())[conf] = previous_energy_profile_[conf];
      }
    }
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
  feasst_serialize(cycles_to_complete_, ostr);
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
  feasst_deserialize(&cycles_to_complete_, istr);
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
  if (config >= static_cast<int>(previous_energy_.size())) {
    previous_energy_.resize(config + 1);
  }
  previous_energy_[config] = current_energy(config);
  (*current_energy_())[config] = energy;
  DEBUG("previous " << feasst_str(previous_energy_));
}

void Criteria::set_current_energy_profile(const std::vector<double>& energy,
    const int config) {
  if (config == static_cast<int>(current_energy_profile_()->size())) {
    current_energy_profile_()->resize(config + 1);
  }
  if (config >= static_cast<int>(previous_energy_profile_.size())) {
    previous_energy_profile_.resize(config + 1);
  }
  previous_energy_profile_[config] = current_energy_profile(config);
  (*current_energy_profile_())[config] = energy;
  DEBUG("previous profile conf" << config << " " << feasst_str(previous_energy_profile_[config]));
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

void Criteria::check_num_cycles_(const int num_trials_per_cycle) {
  *num_attempt_since_last_cycle_() += 1;
  if (*num_attempt_since_last_cycle_() >= num_trials_per_cycle) {
    *num_cycles_() += 1;
    *num_attempt_since_last_cycle_() = 0;
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

int Criteria::num_cycles(const int state) const {
  return const_num_cycles_();
}

void Criteria::initialize(System * system) {
  for (int iconf = 0; iconf < system->num_configurations(); ++iconf) {
    const double en = system->initialize(iconf);
    // HWH set up a Criteria::precompute for this instead.
    set_current_energy(en, iconf);
    set_current_energy_profile(system->stored_energy_profile(iconf), iconf);
  }
  precompute(system);
  update_state(*system, Acceptance());
}

void Criteria::synchronize_(const Criteria& criteria) {
  data_ = criteria.data();
}

}  // namespace feasst
