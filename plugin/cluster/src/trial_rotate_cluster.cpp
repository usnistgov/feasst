#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/trial_rotate_cluster.h"

namespace feasst {

class MapTrialRotateCluster {
 public:
  MapTrialRotateCluster() {
    auto obj = MakeTrialRotateCluster();
    obj->deserialize_map()["TrialRotateCluster"] = obj;
  }
};

static MapTrialRotateCluster mapper_ = MapTrialRotateCluster();

TrialRotateCluster::TrialRotateCluster(argtype * args) : Trial(args) {
  class_name_ = "TrialRotateCluster";
  set_description("TrialRotateCluster");
  add_stage(std::make_shared<SelectCluster>(args),
            std::make_shared<PerturbRotateCOM>(args));
  set(std::make_shared<ComputeMoveCluster>());
}
TrialRotateCluster::TrialRotateCluster(argtype args) : TrialRotateCluster(&args) {
  feasst_check_all_used(args);
}

TrialRotateCluster::TrialRotateCluster(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4538, "mismatch version: " << version);
}

void TrialRotateCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(4538, ostr);
}

}  // namespace feasst
