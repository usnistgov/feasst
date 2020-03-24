#include "steppers/include/num_particles.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapNumParticles {
 public:
  MapNumParticles() {
    auto obj = MakeNumParticles();
    obj->deserialize_map()["NumParticles"] = obj;
  }
};

static MapNumParticles mapper_ = MapNumParticles();

NumParticles::NumParticles(const argtype &args) : Analyze(args) {
  args_.init(args);
  num_particles_.set_block(args_.key("num_blocks").dflt(str(1e5)).integer());
  particle_type_ = args_.key("particle_type").dflt("-1").integer();
  group_ = args_.key("group").dflt("-1").integer();
}

void NumParticles::update(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  ASSERT(particle_type_ == -1 || group_ == -1,
    "both particle type(" << particle_type_ << ") and group(" << group_ <<
    ") cannot be specified at the same time.");
  const Configuration& config = system.configuration();
  DEBUG(particle_type_);
  DEBUG(group_);
  if (particle_type_ == -1) {
    if (group_ == -1) {
      num_particles_.accumulate(config.num_particles(0));
    } else {
      num_particles_.accumulate(config.num_particles(group_));
    }
  } else {
    num_particles_.accumulate(config.num_particles_of_type(particle_type_));
  }
}

std::string NumParticles::write(const Criteria * criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << num_particles_.str() << " ";
  DEBUG(ss.str());
  return ss.str();
}

void NumParticles::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(956, ostr);
  feasst_serialize_fstobj(num_particles_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(group_, ostr);
}

NumParticles::NumParticles(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 956, "version mismatch:" << version);
  feasst_deserialize_fstobj(&num_particles_, istr);
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&group_, istr);
}

}  // namespace feasst
