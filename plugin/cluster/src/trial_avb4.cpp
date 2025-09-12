#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

void gen_avb4_args_(argtype * args, argtype * perturb_args) {
  args->insert({"grand_canonical", "false"});
  args->insert({"second_target", "true"});
  args->insert({"inside", "true"});
  if (used("neighbor_index", *args)) {
    const std::string nindex = str("neighbor_index", args);
    for (argtype * args : {args, perturb_args}) {
      args->insert({"neighbor_index", nindex});
    }
  }
  DEBUG("perturb_args " << str(*perturb_args));
  DEBUG("args " << str(*args));
}

FEASST_MAPPER(TrialAVB4,);

TrialAVB4::TrialAVB4(argtype * args) : Trial(args) {
  class_name_ = "TrialAVB4";
  argtype perturb_args;
  gen_avb4_args_(args, &perturb_args);
  set_description("TrialAVB4");
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbMoveAVB>(perturb_args),
    args
  );
  set(MakeComputeAVB4());
}
TrialAVB4::TrialAVB4(argtype args) : TrialAVB4(&args) {
  feasst_check_all_used(args);
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
