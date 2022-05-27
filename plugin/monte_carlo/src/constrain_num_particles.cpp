#include <cmath>
#include "utils/include/serialize.h"
#include "monte_carlo/include/constrain_num_particles.h"

namespace feasst {

ConstrainNumParticles::ConstrainNumParticles(argtype * args) : Constraint() {
  class_name_ = "ConstrainNumParticles";
  maximum_ = integer("maximum", args, -1);
  minimum_ = integer("minimum", args, 0);
  type_ = integer("type", args, -1);
}
ConstrainNumParticles::ConstrainNumParticles(argtype args) : ConstrainNumParticles(&args) {
  FEASST_CHECK_ALL_USED(args);
}

int ConstrainNumParticles::num_particles(const System& system,
    const Acceptance& acceptance) const {
  const int shift = acceptance.macrostate_shift();
  int num = -1;
  if (type_ != -1) {
    if (acceptance.macrostate_shift_type() == type_ ||
        acceptance.macrostate_shift_type() == -1) {
      num = system.configuration().num_particles_of_type(type_) + shift;
    } else {
      num = system.configuration().num_particles_of_type(type_);
    }
  } else {
    num = system.configuration().num_particles() + shift;
  }
  DEBUG("type: " << type_ << " shift " << shift << " mst " << acceptance.macrostate_shift_type());
  DEBUG("num " << num << " conf.num " << system.configuration().num_particles());
  return num;
}

bool ConstrainNumParticles::is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
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

class MapConstrainNumParticles {
 public:
  MapConstrainNumParticles() {
    auto obj = MakeConstrainNumParticles();
    obj->deserialize_map()["ConstrainNumParticles"] = obj;
  }
};

static MapConstrainNumParticles mapper_ = MapConstrainNumParticles();

ConstrainNumParticles::ConstrainNumParticles(std::istream& istr)
  : Constraint(istr) {
  // ASSERT(class_name_ == "ConstrainNumParticles", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(2492 == version, "mismatch version: " << version);
  feasst_deserialize(&maximum_, istr);
  feasst_deserialize(&minimum_, istr);
  feasst_deserialize(&type_, istr);
}

void ConstrainNumParticles::serialize_constrain_num_particles_(
    std::ostream& ostr) const {
  serialize_constraint_(ostr);
  feasst_serialize_version(2492, ostr);
  feasst_serialize(maximum_, ostr);
  feasst_serialize(minimum_, ostr);
  feasst_serialize(type_, ostr);
}

void ConstrainNumParticles::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_constrain_num_particles_(ostr);
}

}  // namespace feasst
