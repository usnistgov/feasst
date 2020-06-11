#include "utils/include/serialize.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

TrialTranslateCluster::TrialTranslateCluster(
  std::shared_ptr<NeighborCriteria> neighbor_criteria,
  const argtype& args)
  : Trial(args) {
  add_stage(MakeSelectCluster(neighbor_criteria, args),
            MakePerturbTranslate(args));
  set(std::make_shared<ComputeMoveCluster>());
  class_name_ = "TrialTranslateCluster";
}

class MapTrialTranslateCluster {
 public:
  MapTrialTranslateCluster() {
    auto obj = MakeTrialTranslateCluster(MakeNeighborCriteria());
    obj->deserialize_map()["TrialTranslateCluster"] = obj;
  }
};

static MapTrialTranslateCluster mapper_trial_translate_cluster_ = MapTrialTranslateCluster();

std::shared_ptr<Trial> TrialTranslateCluster::create(std::istream& istr) const {
  return std::make_shared<TrialTranslateCluster>(istr);
}

TrialTranslateCluster::TrialTranslateCluster(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialTranslateCluster", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(349 == version, "mismatch version: " << version);
}

void TrialTranslateCluster::serialize_trial_translate_cluster_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(349, ostr);
}

void TrialTranslateCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_translate_cluster_(ostr);
}

TrialRotateCluster::TrialRotateCluster(
  std::shared_ptr<NeighborCriteria> neighbor_criteria,
  const argtype& args)
  : Trial(args) {
  add_stage(MakeSelectCluster(neighbor_criteria, args),
            MakePerturbRotateCOM(args));
  set(std::make_shared<ComputeMoveCluster>());
  class_name_ = "TrialRotateCluster";
}

class MapTrialRotateCluster {
 public:
  MapTrialRotateCluster() {
    auto obj = MakeTrialRotateCluster(MakeNeighborCriteria());
    obj->deserialize_map()["TrialRotateCluster"] = obj;
  }
};

static MapTrialRotateCluster mapper_trial_rotate_cluster_ = MapTrialRotateCluster();

std::shared_ptr<Trial> TrialRotateCluster::create(std::istream& istr) const {
  return std::make_shared<TrialRotateCluster>(istr);
}

TrialRotateCluster::TrialRotateCluster(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialRotateCluster", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(349 == version, "mismatch version: " << version);
}

void TrialRotateCluster::serialize_trial_rotate_cluster_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(349, ostr);
}

void TrialRotateCluster::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_rotate_cluster_(ostr);
}

TrialRigidCluster::TrialRigidCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : TrialFactory() {
  Arguments args_(args);
  args_.dont_check();
  ASSERT(!args_.key("tunable_param").used(),
    "tunable_param args should not be used in this constructor");
  set_weight(args_.key("weight").dflt("1.").dble());
  argtype rot_args = args;
  argtype trans_args = args;
  trans_args.insert({"tunable_param",
    args_.key("translate_param").dflt("0.1").str()});
  rot_args.insert({"tunable_param",
    args_.key("rotate_param").dflt("25").str()});
  add(MakeTrialTranslateCluster(neighbor_criteria, trans_args));
  add(MakeTrialRotateCluster(neighbor_criteria, rot_args));
}

}  // namespace feasst
