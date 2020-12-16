#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

void gen_avb4_args_(const argtype& args, argtype * args_sel, argtype * args_mv) {
  *args_sel = args;
  args_sel->insert({"grand_canonical", "false"});
  args_sel->insert({"second_target", "true"});
  *args_mv = args;
  args_mv->insert({"inside", "true"});
}

std::shared_ptr<Trial> MakeTrialAVB4(const argtype &args) {
  argtype args_sel, args_mv;
  gen_avb4_args_(args, &args_sel, &args_mv);
  auto trial = MakeTrial(args);
  trial->set_description("TrialAVB4");
  trial->add_stage(
    MakeSelectParticleAVB(args_sel),
    MakePerturbMoveAVB(args_mv),
    args
  );
  trial->set(MakeComputeAVB4());
  return trial;
}

}  // namespace feasst
