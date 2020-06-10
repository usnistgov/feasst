#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/trial_avb2_half.h"

namespace feasst {

TrialAVB2Half::TrialAVB2Half(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialAVB2Half";

  argtype args_sel(args);
  args_sel.insert({"grand_canonical", "false"});

  argtype args_mv(args);

  auto compute = MakeComputeAVB2(args);

  Arguments args_(args);
  args_.dont_check();
  if (args_.key("out_to_in").boolean()) {
    args_sel.insert({"inside", "false"});
    args_mv.insert({"inside", "true"});
    set_weight(weight()*compute->probability_out_to_in());
  } else {
    args_sel.insert({"inside", "true"});
    args_mv.insert({"inside", "false"});
    set_weight(weight()*(1. - compute->probability_out_to_in()));
  }

  add_stage(
    MakeSelectParticleAVB(neighbor_criteria, args_sel),
    MakePerturbMoveAVB(neighbor_criteria, args_mv),
    args
  );
  set(compute);
}

class MapTrialAVB2Half {
 public:
  MapTrialAVB2Half() {
    auto obj = MakeTrialAVB2Half(MakeNeighborCriteria(), {{"out_to_in", "true"}});
    obj->deserialize_map()["TrialAVB2Half"] = obj;
  }
};

static MapTrialAVB2Half mapper_ = MapTrialAVB2Half();

std::shared_ptr<Trial> TrialAVB2Half::create(std::istream& istr) const {
  return std::make_shared<TrialAVB2Half>(istr);
}

TrialAVB2Half::TrialAVB2Half(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAVB2Half", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(9274 == version, "mismatch version: " << version);
}

void TrialAVB2Half::serialize_trial_avb2_half_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(9274, ostr);
}

void TrialAVB2Half::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_avb2_half_(ostr);
}

}  // namespace feasst
