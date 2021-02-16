#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

void gen_avb4_args_(argtype * args) {
  args->insert({"grand_canonical", "false"});
  args->insert({"second_target", "true"});
  args->insert({"inside", "true"});
}

std::shared_ptr<Trial> MakeTrialAVB4(argtype args) {
  auto trial = std::make_shared<Trial>(&args);
  gen_avb4_args_(&args);
  trial->set_description("TrialAVB4");
  trial->add_stage(
    std::make_shared<SelectParticleAVB>(&args),
    std::make_shared<PerturbMoveAVB>(&args),
    &args
  );
  trial->set(MakeComputeAVB4());
  check_all_used(args);
  return trial;
}

}  // namespace feasst
