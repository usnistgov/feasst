#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "chain/include/select_particle_pivot.h"
#include "chain/include/perturb_particle_pivot.h"
#include "chain/include/trial_particle_pivot.h"

namespace feasst {

class MapTrialParticlePivot {
 public:
  MapTrialParticlePivot() {
    auto obj = MakeTrialParticlePivot();
    obj->deserialize_map()["TrialParticlePivot"] = obj;
  }
};

static MapTrialParticlePivot mapper_ = MapTrialParticlePivot();

TrialParticlePivot::TrialParticlePivot(argtype * args) :
  TrialMove(std::make_shared<SelectParticlePivot>(args),
            std::make_shared<PerturbParticlePivot>(args),
            args) {
  class_name_ = "TrialParticlePivot";
  set_description("TrialParticlePivot");
}
TrialParticlePivot::TrialParticlePivot(argtype args) : TrialParticlePivot(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialParticlePivot::TrialParticlePivot(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8790, "mismatch version: " << version);
}

void TrialParticlePivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(8790, ostr);
}

}  // namespace feasst
