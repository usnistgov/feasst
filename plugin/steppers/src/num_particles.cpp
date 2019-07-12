#include "steppers/include/num_particles.h"

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

}  // namespace feasst
