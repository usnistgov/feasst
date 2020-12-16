#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_remove.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/trial_transfer_avb.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAddAVB(const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialAddAVB");
  argtype args_sel(args);
  args_sel.insert({"ghost", "true"});
  args_sel.insert({"grand_canonical", "true"});
  trial->add_stage(
    MakeSelectParticleAVB(args_sel),
    MakePerturbAddAVB(args),
    args);
  trial->set(MakeComputeAddAVB());
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemoveAVB(const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialRemoveAVB");
  argtype args_sel(args);
  args_sel.insert({"load_coordinates", "false"});
  args_sel.insert({"grand_canonical", "true"});
  trial->add_stage(
    MakeSelectParticleAVB(args_sel),
    MakePerturbRemove(),
    args);
  trial->set(MakeComputeRemoveAVB());
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransferAVB(const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  factory->add(MakeTrialAddAVB(args));
  factory->add(MakeTrialRemoveAVB(args));
  return factory;
}

}  // namespace feasst
