#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/perturb_add.h"
#include "ewald/include/trial_add_pair.h"

namespace feasst {

TrialAddPair::TrialAddPair(const argtype& args) : Trial(args) {
  // Separate particle_type0 and particle_type1 for respective stages.
  Arguments args_(args);
  args_.dont_check();
  const std::string pt0 = args_.key("particle_type0").remove().str();
  const std::string pt1 = args_.key("particle_type1").remove().str();
  argtype args0 = args_.args();
  argtype args1 = args_.args();
  args0.insert({"particle_type", pt0});
  args1.insert({"particle_type", pt1});

  add_stage(
    std::make_shared<TrialSelectParticle>(args0),
    std::make_shared<PerturbAdd>(args),
    args
  );
  add_stage(
    std::make_shared<TrialSelectParticle>(args1),
    std::make_shared<PerturbAdd>(args),
    args
  );
  set(std::make_shared<TrialComputeAdd>());
  class_name_ = "TrialAddPair";
}

class MapTrialAddPair {
 public:
  MapTrialAddPair() {
    auto obj = MakeTrialAddPair({{"particle_type0", "0"},
                                 {"particle_type1", "1"}});
    obj->deserialize_map()["TrialAddPair"] = obj;
  }
};

static MapTrialAddPair mapper_ = MapTrialAddPair();

std::shared_ptr<Trial> TrialAddPair::create(std::istream& istr) const {
  return std::make_shared<TrialAddPair>(istr);
}

TrialAddPair::TrialAddPair(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAddPair", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(736 == version, "mismatch version: " << version);
}

void TrialAddPair::serialize_trial_add_pair_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(736, ostr);
}

void TrialAddPair::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_add_pair_(ostr);
}

}  // namespace feasst
