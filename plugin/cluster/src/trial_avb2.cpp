#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/trial_avb2.h"

namespace feasst {

TrialAVB2::TrialAVB2(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialAVB2";

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

class MapTrialAVB2 {
 public:
  MapTrialAVB2() {
    auto obj = MakeTrialAVB2(MakeNeighborCriteria(), {{"out_to_in", "true"}});
    obj->deserialize_map()["TrialAVB2"] = obj;
  }
};

static MapTrialAVB2 mapper_ = MapTrialAVB2();

std::shared_ptr<Trial> TrialAVB2::create(std::istream& istr) const {
  return std::make_shared<TrialAVB2>(istr);
}

TrialAVB2::TrialAVB2(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAVB2", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(9274 == version, "mismatch version: " << version);
}

void TrialAVB2::serialize_trial_avb2_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(9274, ostr);
}

void TrialAVB2::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_avb2_(ostr);
}

}  // namespace feasst
