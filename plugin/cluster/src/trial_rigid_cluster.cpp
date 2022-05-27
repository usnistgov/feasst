#include "utils/include/serialize.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

class MapTrialTranslateCluster {
 public:
  MapTrialTranslateCluster() {
    auto obj = MakeTrialTranslateCluster();
    obj->deserialize_map()["TrialTranslateCluster"] = obj;
  }
};

static MapTrialTranslateCluster mapper_trial_translate_cluster_ = MapTrialTranslateCluster();

TrialTranslateCluster::TrialTranslateCluster(argtype * args) : Trial(args) {
  class_name_ = "TrialTranslateCluster";
  set_description("TrialTranslateCluster");
  add_stage(std::make_shared<SelectCluster>(args),
            std::make_shared<PerturbTranslate>(args));
  set(std::make_shared<ComputeMoveCluster>());
}
TrialTranslateCluster::TrialTranslateCluster(argtype args) : TrialTranslateCluster(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialTranslateCluster::TrialTranslateCluster(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4578, "mismatch version: " << version);
}

void TrialTranslateCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(4578, ostr);
}

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
  FEASST_CHECK_ALL_USED(args);
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
