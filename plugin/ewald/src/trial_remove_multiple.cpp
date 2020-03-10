#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "ewald/include/compute_remove_multiple.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_add_multiple.h"

namespace feasst {

TrialRemoveMultiple::TrialRemoveMultiple(const argtype& args) : Trial(args) {
  Arguments args_(args);
  args_.dont_check();
  const std::vector<int> pt = TrialAddMultiple().ptypes(&args_);
  std::vector<argtype> new_args;
  for (int p : pt) {
    argtype nag = args_.args();
    nag.insert({"particle_type", str(p)});
    nag.insert({"load_coordinates", "false"});
    nag.insert({"exclude_perturbed", "true"});
    new_args.push_back(nag);
  }
  for (const argtype& arg : new_args) {
    add_stage(
      MakeTrialSelectParticle(arg),
      MakePerturbRemove(),
      arg);
  }
  set(std::make_shared<ComputeRemoveMultiple>());
  class_name_ = "TrialRemoveMultiple";
}

class MapTrialRemoveMultiple {
 public:
  MapTrialRemoveMultiple() {
    auto obj = MakeTrialRemoveMultiple({{"particle_type0", "0"},
                                    {"particle_type1", "1"}});
    obj->deserialize_map()["TrialRemoveMultiple"] = obj;
  }
};

static MapTrialRemoveMultiple mapper_ = MapTrialRemoveMultiple();

std::shared_ptr<Trial> TrialRemoveMultiple::create(std::istream& istr) const {
  return std::make_shared<TrialRemoveMultiple>(istr);
}

TrialRemoveMultiple::TrialRemoveMultiple(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialRemoveMultiple", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(848 == version, "mismatch version: " << version);
}

void TrialRemoveMultiple::serialize_trial_remove_multiple_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(848, ostr);
}

void TrialRemoveMultiple::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_remove_multiple_(ostr);
}

}  // namespace feasst
