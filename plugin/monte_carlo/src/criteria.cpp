#include <sstream>
#include <iostream>
#include "utils/include/serialize.h"
#include "monte_carlo/include/criteria.h"

namespace feasst {

Criteria::Criteria(const argtype &args) {
  set_trial_state();
  // parse
  args_.init(args);
  if (args_.key("beta").used()) {
    set_beta(args_.dble());
  }
  std::string start("chemical_potential");
  // if only one chemical potential, drop the subscript
  if (args_.key(start).used()) {
    add_chemical_potential(args_.dble());
  } else {
    std::stringstream key;
    int type = 0;
    key << start << type;
    while (args_.key(key.str()).used()) {
      add_chemical_potential(args_.dble());
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  if (args_.key("pH").used()) {
    set_pH(args_.dble());
  }
}

void Criteria::set_beta(const double beta) {
  beta_ = beta;
  beta_initialized_ = true;
}

double Criteria::beta() const {
  ASSERT(beta_initialized_, "beta must be initialized before use");
  return beta_;
}

void Criteria::set_pH(const double pH) {
  pH_ = pH;
  pH_initialized_ = true;
}

double Criteria::pH() const {
  ASSERT(pH_initialized_, "pH must be initialized before use");
  return pH_;
}

double Criteria::chemical_potential(const int particle_type) const {
  ASSERT(particle_type < static_cast<int>(chemical_potentials_.size()),
    "chemical potential of type(" << particle_type <<
    ") must be initalized before use");
  return chemical_potentials_[particle_type];
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
  ss << "beta," << MAX_PRECISION << beta() << std::endl;
  for (int i = 0; i < static_cast<int>(chemical_potentials_.size()); ++i) {
    ss << "mu" << i << "," << MAX_PRECISION << chemical_potentials_[i];
    ss << std::endl;
  }
  return ss.str();
}

void Criteria::set_trial_state(const int state, const int num) {
  trial_state_ = state;
  num_trial_states_ = num;
}

void Criteria::revert(const bool accepted, const double ln_prob) {
  if (accepted) {
    current_energy_ = previous_energy_;
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

std::shared_ptr<Criteria> Criteria::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

bool Criteria::is_equal(const Criteria * criteria) const {
  if (beta_ != criteria->beta_) return false;
  if (current_energy_ != criteria->current_energy_) return false;
  if (trial_state_ != criteria->trial_state_) return false;
  std::stringstream ss1, ss2;
  serialize(ss1);
  criteria->serialize(ss2);
  if (ss1.str() != ss2.str()) {
    INFO(ss1.str());
    INFO(ss2.str());
    return false;
  }
  return true;
}

double Criteria::beta_mu(const int particle_type) const {
  ASSERT(particle_type < static_cast<int>(chemical_potentials_.size()),
    "chemical potential of particle_type: " << particle_type << " not set.");
  return beta()*chemical_potentials_[particle_type];
}

void Criteria::serialize_criteria_(std::ostream& ostr) const {
  feasst_serialize_version(692, ostr);
  feasst_serialize(beta_, ostr);
  feasst_serialize(beta_initialized_, ostr);
  feasst_serialize(pH_, ostr);
  feasst_serialize(pH_initialized_, ostr);
  feasst_serialize(chemical_potentials_, ostr);
  feasst_serialize(current_energy_, ostr);
  feasst_serialize(previous_energy_, ostr);
  feasst_serialize(trial_state_, ostr);
  feasst_serialize(num_trial_states_, ostr);
}

Criteria::Criteria(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 692, "version mismatch: " << version);
  feasst_deserialize(&beta_, istr);
  feasst_deserialize(&beta_initialized_, istr);
  feasst_deserialize(&pH_, istr);
  feasst_deserialize(&pH_initialized_, istr);
  feasst_deserialize(&chemical_potentials_, istr);
  feasst_deserialize(&current_energy_, istr);
  feasst_deserialize(&previous_energy_, istr);
  feasst_deserialize(&trial_state_, istr);
  feasst_deserialize(&num_trial_states_, istr);
}

}  // namespace feasst
