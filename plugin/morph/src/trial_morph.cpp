#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "morph/include/perturb_particle_type.h"
#include "morph/include/compute_morph.h"
#include "morph/include/trial_morph.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialMorph(argtype args) {
  auto trial = MakeTrial(&args);
  trial->set_description("TrialMorph");
  std::string start("particle_type");
  std::stringstream key;
  int type = 0;
  key << start << type;
  //INFO("args: " << str(args));
  argtype stage_args = get_stage_args(&args);
  while (used(key.str(), args)) {
    const std::string type_str = feasst::str(type);
    const std::string ptype = str(key.str(), &args);
    const std::string pmt = str("particle_type_morph" + type_str, &args);
    ASSERT(ptype != pmt, "Attempting to morph particles into the same type: " <<
      ptype);
    auto sel = MakeTrialSelectParticle({{"particle_type", ptype},
                                        {"exclude_perturbed", "true"}});
    argtype tmp_stage_args = stage_args;
    trial->add_stage(sel, MakePerturbParticleType({{"type", pmt}}), &tmp_stage_args);
    ++type;
    ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
    key.str("");
    key << start << type;
  }
  ASSERT(type > 0, "required arguments not used");
  trial->set(MakeComputeMorph());
  check_all_used(args);
  return trial;
}

}  // namespace feasst
