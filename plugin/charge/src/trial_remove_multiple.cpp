#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "charge/include/compute_remove_multiple.h"
#include "charge/include/trial_remove_multiple.h"
#include "charge/include/trial_add_multiple.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialRemoveMultiple(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialRemoveMultiple");
  const std::vector<int> pt = ptypes(&args);
  std::vector<argtype> new_args;
  trial->set(std::make_shared<ComputeRemoveMultiple>(&args));
  for (int p : pt) {
    argtype nag = args;
    nag.insert({"particle_type", str(p)});
    const std::string num_steps = feasst::str("num_steps", &nag, "1");
    nag.insert({"num_steps", num_steps});
    if (num_steps == "1") {
      nag.insert({"load_coordinates", "false"});
    }
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }

  for (argtype arg : new_args) {
    trial->add_stage(
      std::make_shared<TrialSelectParticle>(&arg),
      std::make_shared<PerturbRemove>(),
      &arg);
    check_all_used(arg);
  }
  return trial;
}

}  // namespace feasst
