#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "cluster/include/trial_transfer_avb_divalent.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAddAVBDivalent(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialAddAVBDivalent");
  const std::string particle_type = str("particle_type", &args);
  const std::string particle_type_a = str("particle_type_a", &args);
  const std::string particle_type_b = str("particle_type_b", &args);
  const std::string site_index_a = str("site_index_a", &args, "0");
  const std::string site_index_b = str("site_index_b", &args, "0");
  const std::string neighbor = str("neighbor_index", &args, "0");
  ASSERT(particle_type_a == particle_type_b,
    "hard coded for type_a == type_b in Compute");
  ASSERT(particle_type != particle_type_b,
    "hard coded for type( " << particle_type << ") != type_a(" <<
    particle_type_b << ") in Compute");
  ASSERT(site_index_a == site_index_b,
    "hard coded same");

  // stage0
  argtype sel0_args, parsed_args = args;
  sel0_args.insert({"ghost", "true"});
  sel0_args.insert({"particle_type", particle_type});
  trial->add_stage(
    MakeTrialSelectParticle(sel0_args),
    std::make_shared<PerturbAdd>(&parsed_args),
    &parsed_args);
  check_all_used(parsed_args);

  // stage1
  argtype sel1_args;
  parsed_args = args;
  sel1_args.insert({"ghost", "true"});
  sel1_args.insert({"particle_type", particle_type_a});
  sel1_args.insert({"neighbor_index", neighbor});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(sel1_args),
    std::make_shared<PerturbAddAVB>(&parsed_args),
    &parsed_args);
  check_all_used(parsed_args);

  // stage2
  argtype sel2_args;
  parsed_args = args;
  sel2_args.insert({"ghost", "true"});
  sel2_args.insert({"particle_type", particle_type_b});
  sel2_args.insert({"neighbor_index", neighbor});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(sel2_args),
    std::make_shared<PerturbAddAVB>(&parsed_args),
    &parsed_args);
  check_all_used(parsed_args);

  trial->set(MakeComputeAddAVBDivalent({{"neighbor_index", neighbor}}));
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemoveAVBDivalent(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialRemoveAVBDivalent");
  const std::string particle_type = str("particle_type", &args);
  const std::string particle_type_a = str("particle_type_a", &args);
  const std::string particle_type_b = str("particle_type_b", &args);
  const std::string site_index_a = str("site_index_a", &args, "0");
  const std::string site_index_b = str("site_index_b", &args, "0");
  const std::string neighbor = str("neighbor_index", &args, "0");
  ASSERT(particle_type_a == particle_type_b,
    "hard coded for type_a == type_b in Compute");
  ASSERT(particle_type != particle_type_b,
    "hard coded for type( " << particle_type << ") != type_a(" <<
    particle_type_b << ") in Compute");
  ASSERT(site_index_a == site_index_b,
    "hard coded same");

  // stage0
  argtype sel0_args;
  argtype parsed_args = args;
  sel0_args.insert({"ghost", "false"});
  sel0_args.insert({"particle_type", particle_type});
  trial->add_stage(
    MakeTrialSelectParticle(sel0_args),
    MakePerturbRemove(),
    &parsed_args);
  check_all_used(parsed_args);

  // stage1
  argtype sel1_args;
  parsed_args = args;
  sel1_args.insert({"ghost", "false"});
  sel1_args.insert({"particle_type", particle_type_a});
  sel1_args.insert({"neighbor_index", neighbor});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(sel1_args),
    MakePerturbRemove(),
    &parsed_args);
  check_all_used(parsed_args);

  // stage2
  argtype sel2_args;
  parsed_args = args;
  sel2_args.insert({"ghost", "false"});
  sel2_args.insert({"particle_type", particle_type_b});
  sel2_args.insert({"neighbor_index", neighbor});
  trial->add_stage(
    MakeSelectParticleAVBDivalent(sel2_args),
    MakePerturbRemove(),
    &parsed_args);
  check_all_used(parsed_args);

  trial->set(MakeComputeRemoveAVBDivalent({{"neighbor_index", neighbor}}));
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransferAVBDivalent(
    argtype args) {
  argtype orig_args = args;
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialAddAVBDivalent(orig_args));
  factory->add(MakeTrialRemoveAVBDivalent(orig_args));
  return factory;
}

}  // namespace feasst
