
#include <sstream>
#include "monte_carlo/include/criteria.h"

namespace feasst {

Criteria::Criteria(const argtype &args) {
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
}

void Criteria::set_beta(const double beta) {
  beta_ = beta;
  beta_initialized_ = true;
}

double Criteria::beta() const {
  ASSERT(beta_initialized_, "beta must be initialized before use");
  return beta_;
}

double Criteria::chemical_potential(const int particle_type) const {
  ASSERT(particle_type < static_cast<int>(chemical_potentials_.size()),
    "chemical potential of type(" << particle_type << ") must be initalized before use");
  return chemical_potentials_[particle_type];
}

std::string Criteria::status() const {
  std::stringstream ss;
  ss << current_energy();
  return ss.str();
}

std::string Criteria::write() const {
  std::stringstream ss;
  ss << "beta " << beta() << endl;
  for (double mu : chemical_potentials_) {
    ss << "mu " << mu << endl;
  }
  return ss.str();
}

std::map<std::string, std::shared_ptr<Criteria> >& Criteria::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Criteria> >* ans =
     new std::map<std::string, std::shared_ptr<Criteria> >();
  return *ans;
}

void Criteria::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Criteria> Criteria::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Criteria> Criteria::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void Criteria::serialize_criteria_(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
  feasst_serialize(beta_, ostr);
  feasst_serialize(beta_initialized_, ostr);
  feasst_serialize(chemical_potentials_, ostr);
  feasst_serialize(current_energy_, ostr);
}

Criteria::Criteria(std::istream& istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&beta_, istr);
  feasst_deserialize(&beta_initialized_, istr);
  feasst_deserialize(&chemical_potentials_, istr);
  feasst_deserialize(&current_energy_, istr);
}

}  // namespace feasst
