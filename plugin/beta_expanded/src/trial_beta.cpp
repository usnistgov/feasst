#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "beta_expanded/include/select_nothing.h"
#include "beta_expanded/include/perturb_beta.h"
#include "beta_expanded/include/compute_beta.h"
#include "beta_expanded/include/trial_beta.h"

namespace feasst {

class MapTrialBeta {
 public:
  MapTrialBeta() {
    auto obj = MakeTrialBeta({{"fixed_beta_change", "1"}});
    obj->deserialize_map()["TrialBeta"] = obj;
  }
};

static MapTrialBeta mapper_ = MapTrialBeta();

TrialBeta::TrialBeta(argtype * args) : Trial(args) {
  class_name_ = "TrialBeta";
  set_description("TrialBeta");
  add_stage(
    std::make_shared<SelectNothing>(args),
    std::make_shared<PerturbBeta>(args),
    args);
  set(MakeComputeBeta());
}
TrialBeta::TrialBeta(argtype args) : TrialBeta(&args) {
  feasst_check_all_used(args);
}

TrialBeta::TrialBeta(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3928, "mismatch version: " << version);
}

void TrialBeta::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3928, ostr);
}

}  // namespace feasst
