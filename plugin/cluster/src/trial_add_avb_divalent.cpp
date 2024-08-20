#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "monte_carlo/include/perturb_add.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "cluster/include/trial_add_avb_divalent.h"

namespace feasst {

class MapTrialAddAVBDivalent {
 public:
  MapTrialAddAVBDivalent() {
    auto obj = MakeTrialAddAVBDivalent({{"particle_type", "0"},
                                        {"particle_type_a", "1"},
                                        {"particle_type_b", "1"}});
    obj->deserialize_map()["TrialAddAVBDivalent"] = obj;
  }
};

static MapTrialAddAVBDivalent mapper_ = MapTrialAddAVBDivalent();

TrialAddAVBDivalent::TrialAddAVBDivalent(argtype * args) : Trial(args) {
  class_name_ = "TrialAddAVBDivalent";
  set_description("TrialAddAVBDivalent");
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
  argtype sel0_args, parsed_args = *args;
  sel0_args.insert({"ghost", "true"});
  sel0_args.insert({"particle_type", particle_type});
  parsed_args.insert({"delay_add", "false"});
  auto perturb0 = std::make_shared<PerturbAdd>(&parsed_args);
  ASSERT(!perturb0->delay_add(), "ComputeAddAVBDivalent assumes not delay_add");
  add_stage(
    MakeTrialSelectParticle(sel0_args),
    perturb0,
    &parsed_args);
  feasst_check_all_used(parsed_args);

  // stage1
  argtype sel1_args;
  parsed_args = *args;
  sel1_args.insert({"ghost", "true"});
  sel1_args.insert({"particle_type", particle_type_a});
  sel1_args.insert({"neighbor_index", neighbor});
  parsed_args.insert({"delay_add", "false"});
  auto perturb1 = std::make_shared<PerturbAddAVB>(&parsed_args);
  ASSERT(!perturb1->delay_add(), "ComputeAddAVBDivalent assumes not delay_add");
  add_stage(
    MakeSelectParticleAVBDivalent(sel1_args),
    perturb1,
    &parsed_args);
  feasst_check_all_used(parsed_args);

  // stage2
  argtype sel2_args;
  parsed_args = *args;
  sel2_args.insert({"ghost", "true"});
  sel2_args.insert({"particle_type", particle_type_b});
  sel2_args.insert({"neighbor_index", neighbor});
  parsed_args.insert({"delay_add", "false"});
  auto perturb2 = std::make_shared<PerturbAddAVB>(&parsed_args);
  ASSERT(!perturb2->delay_add(), "ComputeAddAVBDivalent assumes not delay_add");
  add_stage(
    MakeSelectParticleAVBDivalent(sel2_args),
    perturb2,
    &parsed_args);
  feasst_check_all_used(parsed_args);
  set(MakeComputeAddAVBDivalent({{"neighbor_index", neighbor}}));
}
TrialAddAVBDivalent::TrialAddAVBDivalent(argtype args) : TrialAddAVBDivalent(&args) {
  //feasst_check_all_used(args);
}

TrialAddAVBDivalent::TrialAddAVBDivalent(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1345, "mismatch version: " << version);
}

void TrialAddAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(1345, ostr);
}

}  // namespace feasst
