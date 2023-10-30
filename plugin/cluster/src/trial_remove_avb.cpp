#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_remove_avb.h"

namespace feasst {

class MapTrialRemoveAVB {
 public:
  MapTrialRemoveAVB() {
    auto obj = MakeTrialRemoveAVB();
    obj->deserialize_map()["TrialRemoveAVB"] = obj;
  }
};

static MapTrialRemoveAVB mapper_trial_remove_avb_ = MapTrialRemoveAVB();

TrialRemoveAVB::TrialRemoveAVB(argtype * args) : Trial(args) {
  class_name_ = "TrialRemoveAVB";
  set_description("TrialRemoveAVB");
  args->insert({"grand_canonical", "true"});
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbRemove>(),
    args);
  set(MakeComputeRemoveAVB());
}
TrialRemoveAVB::TrialRemoveAVB(argtype args) : TrialRemoveAVB(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialRemoveAVB::TrialRemoveAVB(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3409, "mismatch version: " << version);
}

void TrialRemoveAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3409, ostr);
}

}  // namespace feasst
