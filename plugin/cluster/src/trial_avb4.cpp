#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

TrialAVB4::TrialAVB4(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialAVB4";

  argtype args_sel(args);
  args_sel.insert({"grand_canonical", "false"});
  args_sel.insert({"second_target", "true"});

  argtype args_mv(args);
  args_mv.insert({"inside", "true"});

  add_stage(
    MakeSelectParticleAVB(neighbor_criteria, args_sel),
    MakePerturbMoveAVB(neighbor_criteria, args_mv),
    args
  );
  set(MakeComputeAVB4());
}

class MapTrialAVB4 {
 public:
  MapTrialAVB4() {
    auto obj = MakeTrialAVB4(MakeNeighborCriteria(), {{"out_to_in", "true"}});
    obj->deserialize_map()["TrialAVB4"] = obj;
  }
};

static MapTrialAVB4 mapper_ = MapTrialAVB4();

std::shared_ptr<Trial> TrialAVB4::create(std::istream& istr) const {
  return std::make_shared<TrialAVB4>(istr);
}

TrialAVB4::TrialAVB4(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAVB4", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3957 == version, "mismatch version: " << version);
}

void TrialAVB4::serialize_trial_avb4_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(3957, ostr);
}

void TrialAVB4::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_avb4_(ostr);
}

}  // namespace feasst
