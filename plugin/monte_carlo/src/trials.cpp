#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/perturb_volume.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trials.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialTranslate(const argtype &args) {
  return MakeTrialMove(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbTranslate>(args),
    "TrialTranslate",
    args);
}

std::shared_ptr<Trial> MakeTrialRotate(const argtype &args) {
  return MakeTrialMove(std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbRotate>(args),
    "TrialRotate",
    args);
}

std::shared_ptr<Trial> MakeTrialAdd(const argtype &args) {
  auto trial = std::make_shared<Trial>(args);
  trial->set_description("TrialAdd");
  trial->add_stage(
    MakeTrialSelectParticle(args),
    MakePerturbAdd(args),
    args);
  trial->set(MakeTrialComputeAdd(args));
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemove(const argtype &args) {
  auto trial = std::make_shared<Trial>(args);
  trial->set_description("TrialRemove");
  argtype args2(args);

  // optimization: do not load coordinates if num_steps == 1, by default
  { Arguments tmpargs(args);
    tmpargs.dont_check();
    if (tmpargs.key("num_steps").dflt("1").integer() == 1) {
      if (!tmpargs.key("load_coordinates").used()) {
        args2.insert({"load_coordinates", "false"});
      }
    }
  }
  trial->add_stage(
    std::make_shared<TrialSelectParticle>(args2),
    std::make_shared<PerturbRemove>(),
    args
  );
  trial->set(std::make_shared<TrialComputeRemove>());
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransfer(const argtype &args) {
  auto factory = std::make_shared<TrialFactory>(args);
  factory->add(MakeTrialAdd(args));
  factory->add(MakeTrialRemove(args));
  return factory;
}

std::shared_ptr<Trial> MakeTrialVolume(const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialVolume");
  trial->add_stage(
    MakeTrialSelectParticle(args),
    MakePerturbVolume(args),
    args);
  trial->set(MakeTrialComputeVolume());
  return trial;
}

}  // namespace feasst
