#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

void gen_avb4_args_(argtype * args) {
  args->insert({"grand_canonical", "false"});
  args->insert({"second_target", "true"});
  args->insert({"inside", "true"});
}

class MapTrialAVB4 {
 public:
  MapTrialAVB4() {
    auto obj = MakeTrialAVB4();
    obj->deserialize_map()["TrialAVB4"] = obj;
  }
};

static MapTrialAVB4 mapper_ = MapTrialAVB4();

TrialAVB4::TrialAVB4(argtype * args) : Trial(args) {
  class_name_ = "TrialAVB4";
  gen_avb4_args_(args);
  set_description("TrialAVB4");
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbMoveAVB>(args),
    args
  );
  set(MakeComputeAVB4());
}
TrialAVB4::TrialAVB4(argtype args) : TrialAVB4(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialAVB4::TrialAVB4(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3846, "mismatch version: " << version);
}

void TrialAVB4::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3846, ostr);
}

}  // namespace feasst
