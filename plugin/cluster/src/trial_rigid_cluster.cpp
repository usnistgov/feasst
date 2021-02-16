#include "utils/include/serialize.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/compute_move_cluster.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialTranslateCluster(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialTranslateCluster");
  trial->add_stage(std::make_shared<SelectCluster>(&args),
            std::make_shared<PerturbTranslate>(&args));
  trial->set(std::make_shared<ComputeMoveCluster>());
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialRotateCluster(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialRotateCluster");
  trial->add_stage(std::make_shared<SelectCluster>(&args),
            std::make_shared<PerturbRotateCOM>(&args));
  trial->set(std::make_shared<ComputeMoveCluster>());
  check_all_used(args);
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialRigidCluster(argtype args) {
  ASSERT(!used("tunable_param", args),
    "tunable_param args should not be used in this constructor");
  const std::string tran_param = str("translate_param", &args, "0.1");
  const std::string rot_param = str("rotate_param", &args, "25");
  argtype rot_args = args;
  argtype trans_args = args;
  trans_args.insert({"tunable_param", tran_param});
  rot_args.insert({"tunable_param", rot_param});
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialTranslateCluster(trans_args));
  factory->add(MakeTrialRotateCluster(rot_args));
  return factory;
}

}  // namespace feasst
