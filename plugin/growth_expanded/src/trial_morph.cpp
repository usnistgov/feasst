#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "growth_expanded/include/perturb_particle_type.h"
#include "growth_expanded/include/compute_morph.h"
#include "growth_expanded/include/trial_morph.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialMorph(const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialMorph");
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
    trial->add_stage(sel, MakePerturbParticleType({{"type", pmt}}), args);
    ++type;
    ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
    key.str("");
    key << start << type;
  }
  ASSERT(type > 0, "required arguments not used");
  trial->set(MakeComputeMorph());
  return trial;
}

}  // namespace feasst
