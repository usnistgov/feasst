#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "cluster/include/trial_remove_avb_divalent.h"

namespace feasst {

FEASST_MAPPER(TrialRemoveAVBDivalent, argtype({{"particle_type", "0"},
  {"particle_type_a", "1"}, {"particle_type_b", "1"}}));

TrialRemoveAVBDivalent::TrialRemoveAVBDivalent(argtype * args) : Trial(args) {
  class_name_ = "TrialRemoveAVBDivalent";
  set_description("TrialRemoveAVBDivalent");
  const std::string particle_type = str("particle_type", args);
  const std::string particle_type_a = str("particle_type_a", args);
  const std::string particle_type_b = str("particle_type_b", args);
  const std::string site_index_a = str("site_index_a", args, "0");
  const std::string site_index_b = str("site_index_b", args, "0");
  const std::string neighbor = str("neighbor_index", args, "0");
  ASSERT(particle_type_a == particle_type_b,
    "hard coded for type_a == type_b in Compute");
  ASSERT(particle_type != particle_type_b,
    "hard coded for type( " << particle_type << ") != type_a(" <<
    particle_type_b << ") in Compute");
  ASSERT(site_index_a == site_index_b,
    "hard coded same");

  // stage0
  argtype sel0_args;
  argtype parsed_args = *args;
  sel0_args.insert({"ghost", "false"});
  sel0_args.insert({"particle_type", particle_type});
  add_stage(
    MakeTrialSelectParticle(sel0_args),
    MakePerturbRemove(),
    &parsed_args);
  feasst_check_all_used(parsed_args);

  // stage1
  argtype sel1_args;
  parsed_args = *args;
  sel1_args.insert({"ghost", "false"});
  sel1_args.insert({"particle_type", particle_type_a});
  sel1_args.insert({"neighbor_index", neighbor});
  add_stage(
    MakeSelectParticleAVBDivalent(sel1_args),
    MakePerturbRemove(),
    &parsed_args);
  feasst_check_all_used(parsed_args);

  // stage2
  argtype sel2_args;
  parsed_args = *args;
  sel2_args.insert({"ghost", "false"});
  sel2_args.insert({"particle_type", particle_type_b});
  sel2_args.insert({"neighbor_index", neighbor});
  add_stage(
    MakeSelectParticleAVBDivalent(sel2_args),
    MakePerturbRemove(),
    &parsed_args);
  feasst_check_all_used(parsed_args);
  set(MakeComputeRemoveAVBDivalent({{"neighbor_index", neighbor}}));
}
TrialRemoveAVBDivalent::TrialRemoveAVBDivalent(argtype args) : TrialRemoveAVBDivalent(&args) {
  //feasst_check_all_used(args);
}

TrialRemoveAVBDivalent::TrialRemoveAVBDivalent(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9678, "mismatch version: " << version);
}

void TrialRemoveAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(9678, ostr);
}

}  // namespace feasst
