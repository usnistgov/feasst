#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_add.h"

namespace feasst {

class MapTrialAdd {
 public:
  MapTrialAdd() {
    auto obj = MakeTrialAdd();
    obj->deserialize_map()["TrialAdd"] = obj;
  }
};

static MapTrialAdd mapper_ = MapTrialAdd();

TrialAdd::TrialAdd(argtype * args) : Trial(args) {
  class_name_ = "TrialAdd";
  set_description("TrialAdd");
  auto perturb = std::make_shared<PerturbAdd>(args);
  ASSERT(perturb->delay_add(), "TrialComputeAdd assumes delay_add is true");
  add_stage(std::make_shared<TrialSelectParticle>(args), perturb, args);
  set(std::make_shared<TrialComputeAdd>(args));
}
TrialAdd::TrialAdd(argtype args) : TrialAdd(&args) {
  feasst_check_all_used(args);
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

}  // namespace feasst
