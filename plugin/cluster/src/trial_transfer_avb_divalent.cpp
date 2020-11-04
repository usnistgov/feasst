#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "cluster/include/trial_transfer_avb_divalent.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialAddAVBDivalent");
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
  trial->add_stage(
    MakeTrialSelectParticle(sel0_args),
    MakePerturbAdd(parsed_args),
    parsed_args);

  // stage1
  argtype sel1_args;
  sel1_args.insert({"ghost", "true"});
  sel1_args.insert({"particle_type", particle_type_a});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel1_args),
    MakePerturbAddAVB(neighbor_criteria, parsed_args),
    parsed_args);

  // stage2
  argtype sel2_args;
  sel2_args.insert({"ghost", "true"});
  sel2_args.insert({"particle_type", particle_type_b});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel2_args),
    MakePerturbAddAVB(neighbor_criteria, parsed_args),
    parsed_args);

  trial->set(MakeComputeAddAVBDivalent(neighbor_criteria));
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialRemoveAVBDivalent");
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
  trial->add_stage(
    MakeTrialSelectParticle(sel0_args),
    MakePerturbRemove(),
    parsed_args);

  // stage1
  argtype sel1_args;
  sel1_args.insert({"ghost", "false"});
  sel1_args.insert({"particle_type", particle_type_a});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel1_args),
    MakePerturbRemove(),
    parsed_args);

  // stage2
  argtype sel2_args;
  sel2_args.insert({"ghost", "false"});
  sel2_args.insert({"particle_type", particle_type_b});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(neighbor_criteria, sel2_args),
    MakePerturbRemove(),
    parsed_args);

  trial->set(MakeComputeRemoveAVBDivalent(neighbor_criteria));
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransferAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  factory->add(MakeTrialAddAVBDivalent(neighbor_criteria, args));
  factory->add(MakeTrialRemoveAVBDivalent(neighbor_criteria, args));
  return factory;
}

}  // namespace feasst
