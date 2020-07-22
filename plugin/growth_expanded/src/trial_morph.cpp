#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "growth_expanded/include/perturb_particle_type.h"
#include "growth_expanded/include/compute_morph.h"
#include "growth_expanded/include/trial_morph.h"

namespace feasst {

class MapTrialMorph {
 public:
  MapTrialMorph() {
    auto obj = MakeTrialMorph({{"particle_type0", "0"},
                               {"particle_type_morph0", "1"}});
    obj->deserialize_map()["TrialMorph"] = obj;
  }
};

static MapTrialMorph mapper_trial_morph_ = MapTrialMorph();

TrialMorph::TrialMorph(const argtype& args) : Trial(args) {
  class_name_ = "TrialMorph";
  Arguments args_(args);
  args_.dont_check();
  DEBUG(args_.status());
  std::string start("particle_type");
  std::stringstream key;
  int type = 0;
  key << start << type;
  while (args_.key(key.str()).used()) {
    const std::string type_str = feasst::str(type);
    const std::string ptype = args_.str();
    const std::string pmt = args_.key("particle_type_morph" + type_str).str();
    ASSERT(ptype != pmt, "Attempting to morph particles into the same type: " <<
      ptype);
    auto sel = MakeTrialSelectParticle({{"particle_type", ptype},
                                        {"exclude_perturbed", "true"}});
    add_stage(sel, MakePerturbParticleType({{"type", pmt}}), args);
    ++type;
    ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
    key.str("");
    key << start << type;
  }
  ASSERT(type > 0, "required arguments not used");
  set(MakeComputeMorph());
}

std::shared_ptr<Trial> TrialMorph::create(std::istream& istr) const {
  return std::make_shared<TrialMorph>(istr);
}

TrialMorph::TrialMorph(std::istream& istr) : Trial(istr) {
  ASSERT(class_name_ == "TrialMorph", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(1526 == version, "mismatch version: " << version);
}


void TrialMorph::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(1526, ostr);
}

}  // namespace feasst
