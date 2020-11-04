#include "utils/include/serialize.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialAVB4(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialAVB4");
  argtype args_sel(args);
  args_sel.insert({"grand_canonical", "false"});
  args_sel.insert({"second_target", "true"});
  argtype args_mv(args);
  args_mv.insert({"inside", "true"});
  trial->add_stage(
    MakeSelectParticleAVB(neighbor_criteria, args_sel),
    MakePerturbMoveAVB(neighbor_criteria, args_mv),
    args
  );
  trial->set(MakeComputeAVB4());
  return trial;
}

}  // namespace feasst
