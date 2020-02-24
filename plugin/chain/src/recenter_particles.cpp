#include "chain/include/recenter_particles.h"

namespace feasst {

RecenterParticles::RecenterParticles(const argtype &args)
  : ModifyUpdateOnly(args) {
  group_index_ = args_.key("group_index").dflt("0").integer();
}

class MapRecenterParticles {
 public:
  MapRecenterParticles() {
    RecenterParticles().deserialize_map()["RecenterParticles"] = MakeRecenterParticles();
  }
};

static MapRecenterParticles mapper_ = MapRecenterParticles();

void RecenterParticles::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(448, ostr);
}

RecenterParticles::RecenterParticles(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 448, "version mismatch:" << version);
}

}  // namespace feasst
