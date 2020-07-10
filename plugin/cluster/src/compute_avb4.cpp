#include <cmath>
#include "monte_carlo/include/trial_select.h"
#include "utils/include/serialize.h"
#include "cluster/include/compute_avb4.h"

namespace feasst {

ComputeAVB4::ComputeAVB4() : TrialComputeMove() {
  class_name_ = "ComputeAVB4";
}

class MapComputeAVB4 {
 public:
  MapComputeAVB4() {
    auto obj = MakeComputeAVB4();
    obj->deserialize_map()["ComputeAVB4"] = obj;
  }
};

static MapComputeAVB4 mapper_ = MapComputeAVB4();

void ComputeAVB4::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  TrialComputeMove::perturb_and_acceptance(
    criteria, system, acceptance, stages, random);
  const TrialSelect& select = (*stages)[0]->trial_select();
  acceptance->add_to_ln_metropolis_prob(std::log(select.probability()));
}

std::shared_ptr<TrialCompute> ComputeAVB4::create(std::istream& istr) const {
  return std::make_shared<ComputeAVB4>(istr);
}

ComputeAVB4::ComputeAVB4(std::istream& istr)
  : TrialComputeMove(istr) {
  // ASSERT(class_name_ == "ComputeAVB4", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3674 == version, "mismatch version: " << version);
}

void ComputeAVB4::serialize_compute_avb4_(std::ostream& ostr) const {
  serialize_trial_compute_move_(ostr);
  feasst_serialize_version(3674, ostr);
}

void ComputeAVB4::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_avb4_(ostr);
}

}  // namespace feasst
