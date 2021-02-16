#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_transfer_avb.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAddAVB(argtype args) {
  auto trial = MakeTrialAddAVB(&args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialAddAVB(argtype * args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialAddAVB");
  args->insert({"grand_canonical", "true"});
  trial->add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbAddAVB>(args),
    args);
  trial->set(MakeComputeAddAVB());
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemoveAVB(argtype args) {
  auto trial = MakeTrialRemoveAVB(&args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemoveAVB(argtype * args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialRemoveAVB");
  args->insert({"grand_canonical", "true"});
  trial->add_stage(
    std::make_shared<SelectParticleAVB>(args),
    std::make_shared<PerturbRemove>(),
    args);
  trial->set(MakeComputeRemoveAVB());
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransferAVB(argtype args) {
  argtype orig_args = args;
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialAddAVB(orig_args));
  factory->add(MakeTrialRemoveAVB(orig_args));
  return factory;
}

}  // namespace feasst
