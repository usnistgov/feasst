#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "ewald/include/trial_remove_pair.h"

namespace feasst {

TrialRemovePair::TrialRemovePair(const argtype& args) : Trial(args) {
  // Separate particle_type0 and particle_type1 for respective stages.
  Arguments args_(args);
  args_.dont_check();
  const std::string pt0 = args_.key("particle_type0").remove().str();
  const std::string pt1 = args_.key("particle_type1").remove().str();
  argtype args0 = args_.args();
  argtype args1 = args_.args();
  args0.insert({"particle_type", pt0});
  args0.insert({"load_coordinates", "false"});
  args1.insert({"particle_type", pt1});
  args1.insert({"load_coordinates", "false"});

  add_stage(
    std::make_shared<TrialSelectParticle>(args0),
    std::make_shared<PerturbRemove>(),
    args
  );
  add_stage(
    std::make_shared<TrialSelectParticle>(args1),
    std::make_shared<PerturbRemove>(),
    args
  );
  set(std::make_shared<TrialComputeRemove>());
  class_name_ = "TrialRemovePair";
}

class MapTrialRemovePair {
 public:
  MapTrialRemovePair() {
    auto obj = MakeTrialRemovePair({{"particle_type0", "0"},
                                    {"particle_type1", "1"}});
    obj->deserialize_map()["TrialRemovePair"] = obj;
  }
};

static MapTrialRemovePair mapper_ = MapTrialRemovePair();

std::shared_ptr<Trial> TrialRemovePair::create(std::istream& istr) const {
  return std::make_shared<TrialRemovePair>(istr);
}

TrialRemovePair::TrialRemovePair(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialRemovePair", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(848 == version, "mismatch version: " << version);
}

void TrialRemovePair::serialize_trial_remove_pair_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(848, ostr);
}

void TrialRemovePair::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_remove_pair_(ostr);
}

}  // namespace feasst
