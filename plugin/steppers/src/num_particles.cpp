#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/num_particles.h"

namespace feasst {

FEASST_MAPPER(NumParticles,);

NumParticles::NumParticles(argtype * args) : Analyze(args) {
  particle_type_name_ = str("particle_type", args, "");
  group_ = integer("group", args, -1);
}
NumParticles::NumParticles(argtype args) : NumParticles(&args) {
  feasst_check_all_used(args);
}

std::string NumParticles::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << num_particles().status_header() << std::endl;
  return ss.str();
}

void NumParticles::initialize(MonteCarlo * mc) {
  printer(header(*mc), output_file(mc->criteria()));
  if (particle_type_name_.empty()) {
    particle_type_ = -1;
  } else {
    particle_type_ = configuration(mc->system()).particle_name_to_type(particle_type_name_);
  }
}

void NumParticles::update(const MonteCarlo& mc) {
  ASSERT(particle_type_ == -1 || group_ == -1,
    "both particle type(" << particle_type_ << ") and group(" << group_ <<
    ") cannot be specified at the same time.");
  const Configuration& config = configuration(mc.system());
  DEBUG(particle_type_);
  DEBUG(group_);
  if (particle_type_ == -1) {
    if (group_ == -1) {
      accumulator_->accumulate(config.num_particles(0));
    } else {
      accumulator_->accumulate(config.num_particles(group_));
    }
  } else {
    accumulator_->accumulate(config.num_particles_of_type(particle_type_));
  }
}

std::string NumParticles::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  ss << num_particles().status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void NumParticles::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(957, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(particle_type_name_, ostr);
  feasst_serialize(group_, ostr);
}

NumParticles::NumParticles(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 956 && version <= 957, "version mismatch:" << version);
  feasst_deserialize(&particle_type_, istr);
  if (version >= 957) {
    feasst_deserialize(&particle_type_name_, istr);
  }
  feasst_deserialize(&group_, istr);
}

}  // namespace feasst
