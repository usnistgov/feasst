#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_add_avb.h"

namespace feasst {

class MapTrialAddAVB {
 public:
  MapTrialAddAVB() {
    auto obj = MakeTrialAddAVB();
    obj->deserialize_map()["TrialAddAVB"] = obj;
  }
};

static MapTrialAddAVB mapper_ = MapTrialAddAVB();

// Note that changes here should also be incorported into TrialGrow
TrialAddAVB::TrialAddAVB(argtype * args) : Trial(args) {
  class_name_ = "TrialAddAVB";
  set_description("TrialAddAVB");
  args->insert({"grand_canonical", "true"});
  auto perturb = std::make_shared<PerturbAddAVB>(args);
  ASSERT(perturb->delay_add(), "ComputeAddAVB assumes delay_add is true");
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    perturb,
    args);
  set(MakeComputeAddAVB());
}
TrialAddAVB::TrialAddAVB(argtype args) : TrialAddAVB(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialAddAVB::TrialAddAVB(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3649, "mismatch version: " << version);
}

void TrialAddAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3649, ostr);
}

}  // namespace feasst
