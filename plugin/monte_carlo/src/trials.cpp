#include "utils/include/debug.h"
#include "utils/include/serialize.h"
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

class MapTrialTranslate {
 public:
  MapTrialTranslate() {
    auto obj = MakeTrialTranslate();
    obj->deserialize_map()["TrialTranslate"] = obj;
  }
};

static MapTrialTranslate mapperTranslate_ = MapTrialTranslate();

TrialTranslate::TrialTranslate(argtype * args) :
  TrialMove(std::make_shared<TrialSelectParticle>(args),
            std::make_shared<PerturbTranslate>(args),
            args) {
  class_name_ = "TrialTranslate";
}
TrialTranslate::TrialTranslate(argtype args) : TrialTranslate(&args) {
  check_all_used(args);
}

TrialTranslate::TrialTranslate(std::istream& istr) : TrialMove(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3056, "mismatch version: " << version);
}

void TrialTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(3056, ostr);
}

std::shared_ptr<Trial> MakeTrialRotate(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<TrialSelectParticle>(&args),
    std::make_shared<PerturbRotate>(&args),
    "TrialRotate",
    &args);
  check_all_used(args);
  return trial;
}

class MapTrialAdd {
 public:
  MapTrialAdd() {
    auto obj = MakeTrialAdd();
    obj->deserialize_map()["TrialAdd"] = obj;
  }
};

static MapTrialAdd mapperAdd_ = MapTrialAdd();

TrialAdd::TrialAdd(argtype * args) : Trial(args) {
  class_name_ = "TrialAdd";
  auto perturb = std::make_shared<PerturbAdd>(args);
  ASSERT(perturb->delay_add(), "TrialComputeAdd assumes delay_add is true");
  add_stage(std::make_shared<TrialSelectParticle>(args), perturb, args);
  set(std::make_shared<TrialComputeAdd>(args));
}
TrialAdd::TrialAdd(argtype args) : TrialAdd(&args) {
  check_all_used(args);
}

TrialAdd::TrialAdd(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3056, "mismatch version: " << version);
}

void TrialAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(3056, ostr);
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
