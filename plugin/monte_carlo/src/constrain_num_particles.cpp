#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/constrain_num_particles.h"

namespace feasst {

FEASST_MAPPER(ConstrainNumParticles,);

ConstrainNumParticles::ConstrainNumParticles(argtype * args) : Constraint() {
  class_name_ = "ConstrainNumParticles";
  maximum_ = integer("maximum", args, -1);
  minimum_ = integer("minimum", args, 0);
  type_name_ = str("type", args, "");
  type_ = -2; // initialize if -2
  if (used("configuration_index", *args)) {
    WARN("Deprecated RefPotential::configuration_index->config (see Configuration::name)");
  }
  configuration_index_ = integer("configuration_index", args, 0);
  config_ = str("config", args, "");
}
ConstrainNumParticles::ConstrainNumParticles(argtype args) : ConstrainNumParticles(&args) {
  feasst_check_all_used(args);
}

int ConstrainNumParticles::num_particles(const System& system,
    const Acceptance& acceptance) {
  if (!config_.empty()) {
    if (!config_set_) {
      configuration_index_ = system.configuration_index(config_);
    }
  }
  const int shift = acceptance.macrostate_shift();
  const Configuration& conf = system.configuration(configuration_index_);
  if (type_ == -2) {
    if (type_name_.empty()) {
      type_ = -1;
    } else {
      type_ = conf.particle_name_to_type(type_name_);
    }
  }
  int num = -1;
  if (type_ != -1) {
    if (acceptance.macrostate_shift_type() == type_ ||
        acceptance.macrostate_shift_type() == -1) {
      num = conf.num_particles_of_type(type_) + shift;
    } else {
      num = conf.num_particles_of_type(type_);
    }
  } else {
    num = conf.num_particles() + shift;
  }
  DEBUG("type: " << type_ << " shift " << shift << " mst " << acceptance.macrostate_shift_type());
  DEBUG("num " << num << " conf.num " << conf.num_particles());
  return num;
}

bool ConstrainNumParticles::is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) {
  const int num = num_particles(system, acceptance);
  bool allowed;
  if (num >= minimum_ && ( (maximum_ == -1) || (num <= maximum_) ) ) {
    allowed = true;
  } else {
    allowed = false;
  }
  DEBUG("num: " << num << " allowed: " << allowed << " max " << maximum_ <<
    " min " << minimum_);
  return allowed;
}

ConstrainNumParticles::ConstrainNumParticles(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "ConstrainNumParticles", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 2492 && version <= 2494, "mismatch version: " << version);
  feasst_deserialize(&maximum_, istr);
  feasst_deserialize(&minimum_, istr);
  feasst_deserialize(&type_, istr);
  if (version >= 2493) {
    feasst_deserialize(&type_name_, istr);
    feasst_deserialize(&configuration_index_, istr);
  }
  if (version >= 2494) {
    feasst_deserialize(&config_, istr);
  }
}

void ConstrainNumParticles::serialize_constrain_num_particles_(
    std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(2494, ostr);
  feasst_serialize(maximum_, ostr);
  feasst_serialize(minimum_, ostr);
  feasst_serialize(type_, ostr);
  feasst_serialize(type_name_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(config_, ostr);
}

void ConstrainNumParticles::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_constrain_num_particles_(ostr);
}

}  // namespace feasst
