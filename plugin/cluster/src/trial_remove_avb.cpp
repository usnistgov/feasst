#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_remove_avb.h"

namespace feasst {

TrialRemoveAVB::TrialRemoveAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  argtype args_sel(args);
  args_sel.insert({"load_coordinates", "false"});
  args_sel.insert({"grand_canonical", "true"});
  add_stage(
    MakeSelectParticleAVB(neighbor_criteria, args_sel),
    MakePerturbRemove(),
    args
  );
  set(MakeComputeRemoveAVB());
  class_name_ = "TrialRemoveAVB";
}

class MapTrialRemoveAVB {
 public:
  MapTrialRemoveAVB() {
    auto obj = MakeTrialRemoveAVB(MakeNeighborCriteria());
    obj->deserialize_map()["TrialRemoveAVB"] = obj;
  }
};

static MapTrialRemoveAVB mapper_ = MapTrialRemoveAVB();

std::shared_ptr<Trial> TrialRemoveAVB::create(std::istream& istr) const {
  return std::make_shared<TrialRemoveAVB>(istr);
}

TrialRemoveAVB::TrialRemoveAVB(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialRemoveAVB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(848 == version, "mismatch version: " << version);
}

void TrialRemoveAVB::serialize_trial_remove_avb_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(848, ostr);
}

void TrialRemoveAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_remove_avb_(ostr);
}

}  // namespace feasst
