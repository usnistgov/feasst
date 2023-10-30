#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/trial_translate_cluster.h"
#include "cluster/include/trial_rotate_cluster.h"
#include "cluster/include/trial_rigid_cluster.h"

namespace feasst {

class MapTrialRigidCluster {
 public:
  MapTrialRigidCluster() {
    auto obj = MakeTrialRigidCluster();
    obj->deserialize_map()["TrialRigidCluster"] = obj;
  }
};

static MapTrialRigidCluster mapper_trial_avb2__ = MapTrialRigidCluster();

TrialRigidCluster::TrialRigidCluster(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialRigidCluster";
  ASSERT(!used("tunable_param", *args),
    "tunable_param args should not be used in this constructor");
  const std::string tran_param = str("translate_param", args, "0.1");
  const std::string rot_param = str("rotate_param", args, "25");
  argtype rot_args = *args;
  args->insert({"tunable_param", tran_param});
  rot_args.insert({"tunable_param", rot_param});
  DEBUG("trans_args " << str(*args));
  add(std::make_shared<TrialTranslateCluster>(args));
  DEBUG("rot_args " << str(rot_args));
  add(MakeTrialRotateCluster(rot_args));
  DEBUG("args" << str(*args));
}
TrialRigidCluster::TrialRigidCluster(argtype args) : TrialRigidCluster(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
