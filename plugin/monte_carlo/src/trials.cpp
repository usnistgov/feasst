#include "utils/include/debug.h"
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

std::shared_ptr<Trial> MakeTrialTranslate(argtype args) {
  auto trial = MakeTrialMove(
    std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbTranslate>(&args),
    "TrialTranslate",
    &args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialRotate(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbRotate>(&args),
    "TrialRotate",
    &args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialAdd(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialAdd");
  trial->add_stage(
    std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbAdd>(&args),
    &args);
  trial->set(std::make_shared<TrialComputeAdd>(&args));
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialRemove(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  trial->set_description("TrialRemove");

  // optimization: do not load coordinates if num_steps == 1, by default
  { argtype tmpargs = args;
    if (integer("num_steps", &tmpargs, 1) == 1) {
      if (!used("load_coordinates", tmpargs)) {
        args.insert({"load_coordinates", "false"});
      }
    }
  }
  trial->add_stage(
    std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbRemove>(),
    &args
  );
  trial->set(std::make_shared<TrialComputeRemove>());
  return trial;
}

std::shared_ptr<TrialFactory> MakeTrialTransfer(argtype args) {
  argtype orig_args = args;
  auto factory = std::make_shared<TrialFactory>(&args);
  factory->add(MakeTrialAdd(orig_args));
  factory->add(MakeTrialRemove(orig_args));
  return factory;
}

std::shared_ptr<Trial> MakeTrialVolume(argtype args) {
  auto trial = MakeTrial(&args);
  trial->set_description("TrialVolume");
  trial->add_stage(
    std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbVolume>(&args),
    &args);
  trial->set(MakeTrialComputeVolume());
  return trial;
}

}  // namespace feasst
