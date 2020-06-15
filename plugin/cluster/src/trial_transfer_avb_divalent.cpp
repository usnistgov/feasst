#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "cluster/include/trial_transfer_avb_divalent.h"

namespace feasst {

TrialAddAVBDivalent::TrialAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialAddAVBDivalent";
  Arguments args_(args);
  args_.dont_check();
  const std::string particle_type =
    args_.key("particle_type").remove().str();
  const std::string particle_type_a =
    args_.key("particle_type_a").remove().str();
  const std::string particle_type_b =
    args_.key("particle_type_b").remove().str();
  const std::string site_index_a =
    args_.key("site_index_a").dflt("0").remove().str();
  const std::string site_index_b =
    args_.key("site_index_b").dflt("0").remove().str();
  argtype parsed_args = args_.args();
  ASSERT(particle_type_a == particle_type_b,
    "hard coded for type_a == type_b in Compute");
  ASSERT(particle_type != particle_type_b,
    "hard coded for type( " << particle_type << ") != type_a(" <<
    particle_type_b << ") in Compute");
  ASSERT(site_index_a == site_index_b,
    "hard coded same");

  // stage0
  argtype sel0_args;
  sel0_args.insert({"ghost", "true"});
  sel0_args.insert({"particle_type", particle_type});
  add_stage(
    MakeTrialSelectParticle(sel0_args),
    MakePerturbAdd(parsed_args),
    parsed_args);

  // stage1
  argtype sel1_args;
  sel1_args.insert({"ghost", "true"});
  sel1_args.insert({"particle_type", particle_type_a});
  add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel1_args),
    MakePerturbAddAVB(neighbor_criteria, parsed_args),
    parsed_args);

  // stage2
  argtype sel2_args;
  sel2_args.insert({"ghost", "true"});
  sel2_args.insert({"particle_type", particle_type_b});
  add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel2_args),
    MakePerturbAddAVB(neighbor_criteria, parsed_args),
    parsed_args);


  set(MakeComputeAddAVBDivalent(neighbor_criteria));
}

class MapTrialAddAVBDivalent {
 public:
  MapTrialAddAVBDivalent() {
    auto obj = MakeTrialAddAVBDivalent(MakeNeighborCriteria(), {
      {"particle_type", "0"},
      {"particle_type_a", "1"},
      {"particle_type_b", "1"}});
    obj->deserialize_map()["TrialAddAVBDivalent"] = obj;
  }
};

static MapTrialAddAVBDivalent mapper_ = MapTrialAddAVBDivalent();

std::shared_ptr<Trial> TrialAddAVBDivalent::create(std::istream& istr) const {
  return std::make_shared<TrialAddAVBDivalent>(istr);
}

TrialAddAVBDivalent::TrialAddAVBDivalent(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAddAVBDivalent", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6847 == version, "mismatch version: " << version);
}

void TrialAddAVBDivalent::serialize_trial_add_avb_divalent_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(6847, ostr);
}

void TrialAddAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_add_avb_divalent_(ostr);
}

TrialRemoveAVBDivalent::TrialRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : Trial(args) {
  class_name_ = "TrialRemoveAVBDivalent";
  Arguments args_(args);
  args_.dont_check();
  const std::string particle_type =
    args_.key("particle_type").remove().str();
  const std::string particle_type_a =
    args_.key("particle_type_a").remove().str();
  const std::string particle_type_b =
    args_.key("particle_type_b").remove().str();
  const std::string site_index_a =
    args_.key("site_index_a").dflt("0").remove().str();
  const std::string site_index_b =
    args_.key("site_index_b").dflt("0").remove().str();
  argtype parsed_args = args_.args();
  ASSERT(particle_type_a == particle_type_b,
    "hard coded for type_a == type_b in Compute");
  ASSERT(particle_type != particle_type_b,
    "hard coded for type( " << particle_type << ") != type_a(" <<
    particle_type_b << ") in Compute");
  ASSERT(site_index_a == site_index_b,
    "hard coded same");

  // stage0
  argtype sel0_args;
  sel0_args.insert({"ghost", "false"});
  sel0_args.insert({"particle_type", particle_type});
  add_stage(
    MakeTrialSelectParticle(sel0_args),
    MakePerturbRemove(),
    parsed_args);

  // stage1
  argtype sel1_args;
  sel1_args.insert({"ghost", "false"});
  sel1_args.insert({"particle_type", particle_type_a});
  add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel1_args),
    MakePerturbRemove(),
    parsed_args);

  // stage2
  argtype sel2_args;
  sel2_args.insert({"ghost", "false"});
  sel2_args.insert({"particle_type", particle_type_b});
  add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel2_args),
    MakePerturbRemove(),
    parsed_args);


  set(MakeComputeRemoveAVBDivalent(neighbor_criteria));
}

class MapTrialRemoveAVBDivalent {
 public:
  MapTrialRemoveAVBDivalent() {
    auto obj = MakeTrialRemoveAVBDivalent(MakeNeighborCriteria(), {
      {"particle_type", "0"},
      {"particle_type_a", "1"},
      {"particle_type_b", "1"}});
    obj->deserialize_map()["TrialRemoveAVBDivalent"] = obj;
  }
};

static MapTrialRemoveAVBDivalent mapper2_ = MapTrialRemoveAVBDivalent();

std::shared_ptr<Trial> TrialRemoveAVBDivalent::create(std::istream& istr) const {
  return std::make_shared<TrialRemoveAVBDivalent>(istr);
}

TrialRemoveAVBDivalent::TrialRemoveAVBDivalent(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialRemoveAVBDivalent", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(3058 == version, "mismatch version: " << version);
}

void TrialRemoveAVBDivalent::serialize_trial_remove_avb_divalent_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(3058, ostr);
}

void TrialRemoveAVBDivalent::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_remove_avb_divalent_(ostr);
}

TrialTransferAVBDivalent::TrialTransferAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args) : TrialFactory() {
  Arguments args_(args);
  args_.dont_check();
  set_weight(args_.key("weight").dflt("1.").dble());
  add(MakeTrialAddAVBDivalent(neighbor_criteria, args));
  add(MakeTrialRemoveAVBDivalent(neighbor_criteria, args));
}

}  // namespace feasst
