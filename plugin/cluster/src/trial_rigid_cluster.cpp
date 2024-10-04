#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "cluster/include/trial_translate_cluster.h"
#include "cluster/include/trial_rotate_cluster.h"
#include "cluster/include/trial_rigid_cluster.h"

namespace feasst {

FEASST_MAPPER(TrialRigidCluster,);

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
  feasst_check_all_used(args);
}

}  // namespace feasst
