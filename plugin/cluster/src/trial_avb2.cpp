#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/trial_avb2.h"

namespace feasst {

void gen_avb2_args_(const argtype& args, argtype * args_sel, argtype * args_mv,
    argtype * args_comp) {
  *args_sel = args;
  args_sel->insert({"grand_canonical", "false"});
  //INFO("args_sel: " << str(*args_sel));
  *args_mv = args;
  Arguments args_(args);
  args_.dont_check();
  //*args_comp = args_.remove("neighbor_index", args);
  if (args_.key("out_to_in").boolean()) {
    args_sel->insert({"inside", "false"});
    args_mv->insert({"inside", "true"});
    args_comp->insert({"out_to_in", "true"});
  } else {
    args_sel->insert({"inside", "true"});
    args_mv->insert({"inside", "false"});
    args_comp->insert({"out_to_in", "false"});
  }
}

std::shared_ptr<Trial> MakeTrialAVB2Half(const argtype &args) {
  argtype args_sel, args_mv, args_comp;
  gen_avb2_args_(args, &args_sel, &args_mv, &args_comp);
  auto trial = MakeTrial(args);
  auto compute = MakeComputeAVB2(args_comp);
  Arguments args2_(args);
  args2_.dont_check();
  if (args2_.key("out_to_in").boolean()) {
    trial->set_description("TrialAVB2out_to_in");
    trial->set_weight(trial->weight()*compute->probability_out_to_in());
  } else {
    trial->set_description("TrialAVB2in_to_out");
    trial->set_weight(trial->weight()*(1. - compute->probability_out_to_in()));
  }
  trial->add_stage(
    MakeSelectParticleAVB(args_sel),
    MakePerturbMoveAVB(args_mv),
    args);
  trial->set(compute);
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialAVB2(const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  argtype out2in_args(args);
  out2in_args.insert({"out_to_in", "true"});
  factory->add(MakeTrialAVB2Half(out2in_args));
  argtype in2out_args(args);
  in2out_args.insert({"out_to_in", "false"});
  factory->add(MakeTrialAVB2Half(in2out_args));
  return factory;
}


}  // namespace feasst
