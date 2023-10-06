#include "utils/include/serialize.h"
#include "steppers/include/wrap_particles.h"

namespace feasst {

class MapWrapParticles {
 public:
  MapWrapParticles() {
    WrapParticles().deserialize_map()["WrapParticles"] = MakeWrapParticles();
  }
};

static MapWrapParticles mapper_energy_check_ = MapWrapParticles();

WrapParticles::WrapParticles(argtype * args) : ModifyUpdateOnly(args) {
}
WrapParticles::WrapParticles(argtype args) : WrapParticles(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void WrapParticles::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  Configuration * config = system->get_configuration();
  for (int particle_index : config->group_select(0).particle_indices()) {
    config->wrap_particle(particle_index);
  }
}

void WrapParticles::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(3609, ostr);
}

WrapParticles::WrapParticles(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3609, "version mismatch: " << version);
}

}  // namespace feasst
