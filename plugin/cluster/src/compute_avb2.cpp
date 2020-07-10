#include <cmath>
#include "monte_carlo/include/trial_select.h"
#include "utils/include/serialize.h"
#include "cluster/include/compute_avb2.h"

namespace feasst {

ComputeAVB2::ComputeAVB2(const argtype& args)
  : TrialComputeMove() {
  class_name_ = "ComputeAVB2";
  Arguments args_(args);
  p_bias_ = args_.key("probability_out_to_in").dflt("0.5").dble();
  ASSERT((p_bias_ > 0.) && (p_bias_ < 1.), "probability_out_to_in: " << p_bias_
    << "must be >0 and <1");
  out_to_in_ = args_.key("out_to_in").boolean();
}

class MapComputeAVB2 {
 public:
  MapComputeAVB2() {
    auto obj = MakeComputeAVB2({{"out_to_in", "true"}});
    obj->deserialize_map()["ComputeAVB2"] = obj;
  }
};

static MapComputeAVB2 mapper_ = MapComputeAVB2();

void ComputeAVB2::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  TrialComputeMove::perturb_and_acceptance(
    criteria, system, acceptance, stages, random);
  const TrialSelect& select = (*stages)[0]->trial_select();
  double factor = p_bias_/(1 - p_bias_);
  if (out_to_in_) {
    factor = (1 - p_bias_)/p_bias_;
  }
  acceptance->add_to_ln_metropolis_prob(std::log(factor*select.probability()));
}

std::shared_ptr<TrialCompute> ComputeAVB2::create(std::istream& istr) const {
  return std::make_shared<ComputeAVB2>(istr);
}

ComputeAVB2::ComputeAVB2(std::istream& istr)
  : TrialComputeMove(istr) {
  // ASSERT(class_name_ == "ComputeAVB2", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(4203 == version, "mismatch version: " << version);
  feasst_deserialize(&p_bias_, istr);
  feasst_deserialize(&out_to_in_, istr);
}

void ComputeAVB2::serialize_compute_avb2_(std::ostream& ostr) const {
  serialize_trial_compute_move_(ostr);
  feasst_serialize_version(4203, ostr);
  feasst_serialize(p_bias_, ostr);
  feasst_serialize(out_to_in_, ostr);
}

void ComputeAVB2::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_compute_avb2_(ostr);
}

}  // namespace feasst
