#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_all.h"
#include "model_expanded/include/perturb_model.h"
#include "model_expanded/include/compute_model.h"
#include "model_expanded/include/trial_model.h"

namespace feasst {

class MapTrialModel {
 public:
  MapTrialModel() {
    auto obj = MakeTrialModel();
    obj->deserialize_map()["TrialModel"] = obj;
  }
};

static MapTrialModel mapper_ = MapTrialModel();

TrialModel::TrialModel(argtype * args) : Trial(args) {
  class_name_ = "TrialModel";
  set_description("TrialModel");
  add_stage(
    std::make_shared<TrialSelectAll>(args),
    std::make_shared<PerturbModel>(args),
    args);
  set(MakeComputeModel());
}
TrialModel::TrialModel(argtype args) : TrialModel(&args) {
  feasst_check_all_used(args);
}

TrialModel::TrialModel(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6792, "mismatch version: " << version);
}

void TrialModel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(6792, ostr);
}

}  // namespace feasst
