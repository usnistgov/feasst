
#include <sstream>
#include "core/include/criteria.h"

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
  ss << running_energy();
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

}  // namespace feasst
