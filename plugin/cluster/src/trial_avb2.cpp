#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/trial_avb2.h"

namespace feasst {

void gen_avb2_args_(const bool out_to_in, argtype * args,
    argtype * perturb_args) {
  args->insert({"grand_canonical", "false"});
  if (out_to_in) {
    args->insert({"inside", "false"}); // sel
    perturb_args->insert({"inside", "true"});  // mv
    args->insert({"out_to_in", "true"});
  } else {
    args->insert({"inside", "true"});  // sel
    perturb_args->insert({"inside", "false"}); // mv
    args->insert({"out_to_in", "false"});
  }
  if (used("neighbor_index", *args)) {
    perturb_args->insert({"neighbor_index", str("neighbor_index", args)});
  }
}

std::shared_ptr<Trial> MakeTrialAVB2Half(argtype args) {
  // read out_to_in and then put it back in args
  const bool out_to_in = boolean("out_to_in", &args);
  if (out_to_in) {
    args.insert({"out_to_in", "true"});
  } else {
    args.insert({"out_to_in", "false"});
  }

  argtype perturb_args;
  gen_avb2_args_(out_to_in, &args, &perturb_args);
  auto trial = MakeTrial(&args);
  auto compute = std::make_shared<ComputeAVB2>(&args);
  if (out_to_in) {
    trial->set_description("TrialAVB2out_to_in");
    trial->set_weight(trial->weight()*compute->probability_out_to_in());
  } else {
    trial->set_description("TrialAVB2in_to_out");
    trial->set_weight(trial->weight()*(1. - compute->probability_out_to_in()));
  }
  trial->add_stage(
    std::make_shared<SelectParticleAVB>(&args),
    std::make_shared<PerturbMoveAVB>(&perturb_args),
    &args);
  trial->set(compute);
  check_all_used(args);
  check_all_used(perturb_args);
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialAVB2(argtype args) {
  argtype out2in_args(args);
  out2in_args.insert({"out_to_in", "true"});
  argtype in2out_args(args);
  in2out_args.insert({"out_to_in", "false"});
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialAVB2Half(out2in_args));
  factory->add(MakeTrialAVB2Half(in2out_args));
  return factory;
}


}  // namespace feasst
