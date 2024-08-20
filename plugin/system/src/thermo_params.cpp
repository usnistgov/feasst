#include "utils/include/utils.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "system/include/thermo_params.h"

namespace feasst {

ThermoParams::ThermoParams(argtype args) {
  if (used("beta", args)) {
    set_beta(dble("beta", &args));
  }
  std::string start("chemical_potential");
  // if only one chemical potential, drop the subscript
  if (used(start, args)) {
    add_chemical_potential(dble(start, &args));
  } else {
    std::stringstream key;
    int type = 0;
    key << start << type;
    while (used(key.str(), args)) {
      add_chemical_potential(dble(key.str(), &args));
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  if (used("pH", args)) set_pH(dble("pH", &args));
  if (used("pressure", args)) set_pressure(dble("pressure", &args));
  feasst_check_all_used(args);
}

void ThermoParams::set_beta(const double beta) {
  beta_ = beta;
  beta_initialized_ = true;
}

double ThermoParams::beta() const {
  ASSERT(beta_initialized_, "beta must be initialized before use");
  return beta_;
}

void ThermoParams::set_pH(const double pH) {
  pH_ = pH;
  pH_initialized_ = true;
}

double ThermoParams::pH() const {
  ASSERT(pH_initialized_, "pH must be initialized before use");
  return pH_;
}

double ThermoParams::chemical_potential(const int particle_type) const {
  ASSERT(particle_type < static_cast<int>(chemical_potentials_.size()),
    "chemical potential of type(" << particle_type <<
    ") must be initalized before use");
  return chemical_potentials_[particle_type];
}

void ThermoParams::set_pressure(const double pressure) {
  pressure_ = pressure;
  pressure_initialized_ = true;
}

double ThermoParams::pressure() const {
  ASSERT(pressure_initialized_, "pressure must be initialized before use");
  return pressure_;
}

std::string ThermoParams::str() const {
  std::stringstream ss;
  if (beta_initialized_) ss << "beta," << MAX_PRECISION << beta() << std::endl;
  if (pressure_initialized_) {
    ss << "pressure," << MAX_PRECISION << pressure() << std::endl;
  }
  if (pH_initialized_) ss << "pH," << MAX_PRECISION << pH() << std::endl;
  for (int i = 0; i < static_cast<int>(chemical_potentials_.size()); ++i) {
    ss << "mu" << i << "," << MAX_PRECISION << chemical_potentials_[i];
    ss << std::endl;
  }
  return ss.str();
}

bool ThermoParams::is_equal(const ThermoParams& thermo_params) const {
  if (beta_ != thermo_params.beta_) return false;
  if (!feasst::is_equal(chemical_potentials_, thermo_params.chemical_potentials_)) {
    return false;
  }
  if (pressure_ != thermo_params.pressure_) return false;
  if (pH_ != thermo_params.pH_) return false;
  return true;
}

double ThermoParams::beta_mu(const int particle_type) const {
  ASSERT(particle_type < static_cast<int>(chemical_potentials_.size()),
    "chemical potential of particle_type: " << particle_type << " not set.");
  return beta()*chemical_potentials_[particle_type];
}

void ThermoParams::serialize(std::ostream& ostr) const {
  feasst_serialize_version(9486, ostr);
  feasst_serialize(beta_, ostr);
  feasst_serialize(beta_initialized_, ostr);
  feasst_serialize(pH_, ostr);
  feasst_serialize(pH_initialized_, ostr);
  feasst_serialize(chemical_potentials_, ostr);
  feasst_serialize(pressure_initialized_, ostr);
  feasst_serialize(pressure_, ostr);
}

ThermoParams::ThermoParams(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9486, "version mismatch: " << version);
  feasst_deserialize(&beta_, istr);
  feasst_deserialize(&beta_initialized_, istr);
  feasst_deserialize(&pH_, istr);
  feasst_deserialize(&pH_initialized_, istr);
  feasst_deserialize(&chemical_potentials_, istr);
  feasst_deserialize(&pressure_initialized_, istr);
  feasst_deserialize(&pressure_, istr);
}

}  // namespace feasst
