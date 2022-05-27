#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "charge/include/compute_remove_multiple.h"
#include "charge/include/trial_remove_multiple.h"
#include "charge/include/trial_add_multiple.h"

namespace feasst {

class MapTrialRemoveMultiple {
 public:
  MapTrialRemoveMultiple() {
    auto obj = MakeTrialRemoveMultiple();
    obj->deserialize_map()["TrialRemoveMultiple"] = obj;
  }
};

static MapTrialRemoveMultiple mapper_ = MapTrialRemoveMultiple();

TrialRemoveMultiple::TrialRemoveMultiple(argtype * args) : Trial(args) {
  class_name_ = "TrialRemoveMultiple";
  set_description("TrialRemoveMultiple");
  const std::vector<int> pt = ptypes(args);
  std::vector<argtype> new_args;
  set(std::make_shared<ComputeRemoveMultiple>(args));
  const std::string num_steps = feasst::str("num_steps", args, "1");
  const std::string reference_index = feasst::str("reference_index", args, "-1");
  for (int p : pt) {
    argtype nag = *args;
    nag.insert({"particle_type", str(p)});
    nag.insert({"num_steps", num_steps});
    if (num_steps == "1") {
      nag.insert({"load_coordinates", "false"});
    }
    nag.insert({"reference_index", reference_index});
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }

  for (argtype arg : new_args) {
    add_stage(
      std::make_shared<TrialSelectParticle>(&arg),
      std::make_shared<PerturbRemove>(),
      &arg);
    FEASST_CHECK_ALL_USED(arg);
  }
}
TrialRemoveMultiple::TrialRemoveMultiple(argtype args) : TrialRemoveMultiple(&args) {
  //FEASST_CHECK_ALL_USED(args);
}

TrialRemoveMultiple::TrialRemoveMultiple(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9586, "mismatch version: " << version);
}

void TrialRemoveMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(9586, ostr);
}

}  // namespace feasst
