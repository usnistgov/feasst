#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "ewald/include/compute_remove_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_add_multiple.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialRemoveMultiple(
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialRemoveMultiple");
  Arguments args_(args);
  args_.dont_check();
  const std::vector<int> pt = ptypes(&args_);
  std::vector<argtype> new_args;
  for (int p : pt) {
    argtype nag = args_.args();
    nag.insert({"particle_type", str(p)});
    if (args_.key("num_steps").dflt("1").integer() == 1) {
      nag.insert({"load_coordinates", "false"});
    }
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }
  for (const argtype& arg : new_args) {
    trial->add_stage(
      MakeTrialSelectParticle(arg),
      MakePerturbRemove(),
      arg);
  }
  trial->set(std::make_shared<ComputeRemoveMultiple>(args));
  return trial;
}

}  // namespace feasst
