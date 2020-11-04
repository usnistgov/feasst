#include "utils/include/serialize.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialTranslateCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialTranslateCluster");
  trial->add_stage(MakeSelectCluster(neighbor_criteria, args),
            MakePerturbTranslate(args));
  trial->set(std::make_shared<ComputeMoveCluster>());
  return trial;
}

std::shared_ptr<Trial> MakeTrialRotateCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialRotateCluster");
  trial->add_stage(MakeSelectCluster(neighbor_criteria, args),
            MakePerturbRotateCOM(args));
  trial->set(std::make_shared<ComputeMoveCluster>());
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialRigidCluster(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  Arguments args_(args);
  args_.dont_check();
  ASSERT(!args_.key("tunable_param").used(),
    "tunable_param args should not be used in this constructor");
  argtype rot_args = args;
  argtype trans_args = args;
  trans_args.insert({"tunable_param",
    args_.key("translate_param").dflt("0.1").str()});
  rot_args.insert({"tunable_param",
    args_.key("rotate_param").dflt("25").str()});
  factory->add(MakeTrialTranslateCluster(neighbor_criteria, trans_args));
  factory->add(MakeTrialRotateCluster(neighbor_criteria, rot_args));
  return factory;
}

}  // namespace feasst
