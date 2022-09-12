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

class MapTrialAVB2Half {
 public:
  MapTrialAVB2Half() {
    auto obj = MakeTrialAVB2Half({{"out_to_in", "true"}});
    obj->deserialize_map()["TrialAVB2Half"] = obj;
  }
};

static MapTrialAVB2Half mapper_ = MapTrialAVB2Half();

TrialAVB2Half::TrialAVB2Half(argtype * args) : Trial(args) {
  class_name_ = "TrialAVB2Half";
  // read out_to_in and then put it back in args
  const bool out_to_in = boolean("out_to_in", args);
  if (out_to_in) {
    args->insert({"out_to_in", "true"});
  } else {
    args->insert({"out_to_in", "false"});
  }

  argtype perturb_args;
  gen_avb2_args_(out_to_in, args, &perturb_args);
  auto compute = std::make_shared<ComputeAVB2>(args);
  if (out_to_in) {
    set_description("TrialAVB2out_to_in");
    set_weight(weight()*compute->probability_out_to_in());
  } else {
    set_description("TrialAVB2in_to_out");
    set_weight(weight()*(1. - compute->probability_out_to_in()));
  }
  add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbMoveAVB>(&perturb_args),
    args);
  set(compute);
//  FEASST_CHECK_ALL_USED(args);
  FEASST_CHECK_ALL_USED(perturb_args);
}
TrialAVB2Half::TrialAVB2Half(argtype args) : TrialAVB2Half(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialAVB2Half::TrialAVB2Half(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1634, "mismatch version: " << version);
}

void TrialAVB2Half::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(1634, ostr);
}

class MapTrialAVB2 {
 public:
  MapTrialAVB2() {
    auto obj = MakeTrialAVB2();
    obj->deserialize_map()["TrialAVB2"] = obj;
  }
};

static MapTrialAVB2 mapper_trial_avb2__ = MapTrialAVB2();

TrialAVB2::TrialAVB2(argtype * args) : TrialFactoryNamed() {
  class_name_ = "TrialAVB2";
  argtype out2in_args(*args);
  out2in_args.insert({"out_to_in", "true"});
  argtype * in2out_args = args;
  in2out_args->insert({"out_to_in", "false"});
  auto trial_out2in = MakeTrialAVB2Half(out2in_args);
  trial_out2in->set_weight(trial_out2in->weight()/2.);
  add(trial_out2in);
  auto trial_in2out = std::make_shared<TrialAVB2Half>(in2out_args);
  trial_in2out->set_weight(trial_in2out->weight()/2.);
  add(trial_in2out);
}
TrialAVB2::TrialAVB2(argtype args) : TrialAVB2(&args) {
  FEASST_CHECK_ALL_USED(args);
}

}  // namespace feasst
