#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/trial_add_avb.h"

namespace feasst {

TrialAddAVB::TrialAddAVB(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  argtype args_sel(args);
  args_sel.insert({"ghost", "true"});
  args_sel.insert({"grand_canonical", "true"});
  add_stage(
    MakeSelectParticleAVB(neighbor_criteria, args_sel),
    MakePerturbAddAVB(neighbor_criteria, args),
    args
  );
  set(MakeComputeAddAVB());
  class_name_ = "TrialAddAVB";
}

class MapTrialAddAVB {
 public:
  MapTrialAddAVB() {
    auto obj = MakeTrialAddAVB(MakeNeighborCriteria());
    obj->deserialize_map()["TrialAddAVB"] = obj;
  }
};

static MapTrialAddAVB mapper_ = MapTrialAddAVB();

std::shared_ptr<Trial> TrialAddAVB::create(std::istream& istr) const {
  return std::make_shared<TrialAddAVB>(istr);
}

TrialAddAVB::TrialAddAVB(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAddAVB", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3055 == version, "mismatch version: " << version);
}

void TrialAddAVB::serialize_trial_add_avb_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(3055, ostr);
}

void TrialAddAVB::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_add_avb_(ostr);
}

}  // namespace feasst
