#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/trial_avb2.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAVB2Half(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto trial = MakeTrial(args);
  argtype args_sel(args);
  args_sel.insert({"grand_canonical", "false"});
  argtype args_mv(args);
  auto compute = MakeComputeAVB2(args);
  Arguments args_(args);
  args_.dont_check();
  if (args_.key("out_to_in").boolean()) {
    trial->set_description("TrialAVB2out_to_in");
    args_sel.insert({"inside", "false"});
    args_mv.insert({"inside", "true"});
    trial->set_weight(trial->weight()*compute->probability_out_to_in());
  } else {
    trial->set_description("TrialAVB2in_to_out");
    args_sel.insert({"inside", "true"});
    args_mv.insert({"inside", "false"});
    trial->set_weight(trial->weight()*(1. - compute->probability_out_to_in()));
  }

  trial->add_stage(
    MakeSelectParticleAVB(neighbor_criteria, args_sel),
    MakePerturbMoveAVB(neighbor_criteria, args_mv),
    args
  );
  trial->set(compute);
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialAVB2(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  argtype out2in_args(args);
  out2in_args.insert({"out_to_in", "true"});
  factory->add(MakeTrialAVB2Half(neighbor_criteria, out2in_args));
  argtype in2out_args(args);
  in2out_args.insert({"out_to_in", "false"});
  factory->add(MakeTrialAVB2Half(neighbor_criteria, in2out_args));
  return factory;
}


}  // namespace feasst
