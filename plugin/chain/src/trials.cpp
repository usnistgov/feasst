#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/select_segment.h"
#include "chain/include/perturb_crankshaft.h"
#include "chain/include/select_reptate.h"
#include "chain/include/perturb_reptate.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/trials.h"
#include "chain/include/select_branch.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialPivot(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<SelectEndSegment>(&args),
    std::make_shared<PerturbPivot>(&args),
    "TrialPivot",
    &args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialCrankshaft(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<SelectSegment>(&args),
    std::make_shared<PerturbCrankshaft>(&args),
    "TrialCrankshaft",
    &args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialReptate(argtype args) {
  auto trial = MakeTrialMove(std::make_shared<SelectReptate>(&args),
    std::make_shared<PerturbReptate>(&args),
    "TrialReptate",
    &args);
  check_all_used(args);
  return trial;
}

std::shared_ptr<Trial> MakeTrialSwapSites(argtype args) {
  auto trial = MakeTrial(&args);
  trial->set_description("TrialSwapSites");
  trial->set(MakeTrialComputeMove());
  const int site_type1 = integer("site_type1", &args);
  const int site_type2 = integer("site_type2", &args);
  ASSERT(site_type1 != site_type2, "site types should not match: " <<
    site_type1 << " " << site_type2);
  const std::string part_type = str("particle_type", &args);
  argtype stage_args = args;
  trial->add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type1)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type2)}}),
    &stage_args);
  check_all_used(stage_args);
  trial->add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type2)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type1)}}),
    &args
  );
  check_all_used(args);
  return trial;
}

}  // namespace feasst
