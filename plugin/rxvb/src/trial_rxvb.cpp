#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "gibbs/include/perturb_particle_type.h"
#include "rxvb/include/select_anchor.h"
#include "rxvb/include/compute_rxvb.h"
#include "rxvb/include/trial_rxvb.h"

namespace feasst {

void gen_rxnavb_args_(const bool avb, argtype * args,
    argtype * perturb_args) {
//  args->insert({"grand_canonical", "false"});
//  if (avb) {
//    args->insert({"inside", "false"}); // sel
//    //perturb_args->insert({"inside", "true"});  // mv
//    args->insert({"avb", "true"});
//  } else {
//    args->insert({"inside", "true"});  // sel
//    //perturb_args->insert({"inside", "false"}); // mv
//    args->insert({"avb", "false"});
//  }
  if (used("neighbor_index", *args)) {
    //perturb_args->insert({"neighbor_index", str("neighbor_index", args)});
  }
}

FEASST_MAPPER(TrialRxVBHalf, argtype({{"avb", "true"}, {"target_particle_type", "0"}, {"target_particle_type_morph", "0"}, {"particle_type_morph", "0"}, {"particle_type", "0"}}));

TrialRxVBHalf::TrialRxVBHalf(argtype * args) : Trial(args) {
  class_name_ = "TrialRxVBHalf";
  // read avb and then put it back in args
  const bool avb = boolean("avb", args);
  if (avb) {
    args->insert({"avb", "true"});
  } else {
    args->insert({"avb", "false"});
  }

//  argtype perturb_args;
//  gen_rxnavb_args_(avb, args, &perturb_args);
  DEBUG("args " << str(*args));
  const std::string target_particle_type_morph = str("target_particle_type_morph", args);
  const std::string particle_type_morph = str("particle_type_morph", args);
  const std::string target_particle_type = str("target_particle_type", args);
  const std::string particle_type = str("particle_type", args);
  const int reference_index = integer("reference_index", args, -1);
  DEBUG("args " << str(*args));
  args->insert({"target_particle_type_morph", target_particle_type_morph});
  args->insert({"particle_type_morph", particle_type_morph});
  args->insert({"target_particle_type", target_particle_type});
  args->insert({"particle_type", particle_type});
  args->insert({"reference_index", str(reference_index)});
  DEBUG("args " << str(*args));
  auto compute = std::make_shared<ComputeRXVB>(args);
  DEBUG("args " << str(*args));
  DEBUG("avb " << avb);
  //ASSERT(args->size() == 0, "error");
  args->insert({"grand_canonical", "false"});
  args->insert({"rxnavb", "true"});
  if (avb) {
    set_description("TrialRxVBforward");
    set_weight(weight()*compute->probability_avb());
    args->insert({"target_particle_type", target_particle_type});
    args->insert({"particle_type", particle_type});
    add_stage(
      std::make_shared<SelectParticleAVB>(args),
      std::make_shared<PerturbParticleType>(argtype({{"type", particle_type_morph}})),
      args);
    args->insert({"reference_index", str(reference_index)});
    add_stage(
      std::make_shared<SelectAnchor>(args),
      std::make_shared<PerturbParticleType>(argtype({{"type", target_particle_type_morph}})),
      args);
  } else {
    set_description("TrialRxVBreverse");
    set_weight(weight()*(1. - compute->probability_avb()));
    args->insert({"target_particle_type", target_particle_type_morph});
    args->insert({"particle_type", particle_type_morph});
    add_stage(
      std::make_shared<SelectParticleAVB>(args),
      std::make_shared<PerturbParticleType>(argtype({{"type", particle_type}})),
      args);
    //args->insert({"particle_type", str(target_particle_type_morph)});
    args->insert({"reference_index", str(reference_index)});
    add_stage(
      std::make_shared<SelectAnchor>(args),
      std::make_shared<PerturbParticleType>(argtype({{"type", target_particle_type}})),
      args);
  }
  set(compute);
//  feasst_check_all_used(perturb_args);
}
TrialRxVBHalf::TrialRxVBHalf(argtype args) : TrialRxVBHalf(&args) {
  feasst_check_all_used(args);
}

TrialRxVBHalf::TrialRxVBHalf(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2967, "mismatch version: " << version);
}

void TrialRxVBHalf::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(2967, ostr);
}

FEASST_MAPPER(TrialRxVB, argtype({{"target_particle_type", "0"}, {"target_particle_type_morph", "0"}, {"particle_type_morph", "0"}, {"particle_type", "0"}}));

TrialRxVB::TrialRxVB(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialRxVB";
  argtype out2in_args(*args);
  out2in_args.insert({"avb", "true"});
  argtype * in2out_args = args;
  in2out_args->insert({"avb", "false"});
  auto trial_out2in = std::make_shared<TrialRxVBHalf>(out2in_args);
  trial_out2in->set_weight(trial_out2in->weight()/2.);
  add(trial_out2in);
  auto trial_in2out = std::make_shared<TrialRxVBHalf>(in2out_args);
  trial_in2out->set_weight(trial_in2out->weight()/2.);
  add(trial_in2out);
}
TrialRxVB::TrialRxVB(argtype args) : TrialRxVB(&args) {
  feasst_check_all_used(args);
}

}  // namespace feasst
