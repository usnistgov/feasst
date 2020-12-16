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

std::shared_ptr<Trial> MakeTrialPivot(const argtype &args) {
  return MakeTrialMove(std::make_shared<SelectEndSegment>(args),
    std::make_shared<PerturbPivot>(args),
    "TrialPivot",
    args);
}

std::shared_ptr<Trial> MakeTrialCrankshaft(const argtype &args) {
  return MakeTrialMove(std::make_shared<SelectSegment>(args),
    std::make_shared<PerturbCrankshaft>(args),
    "TrialCrankshaft",
    args);
}

std::shared_ptr<Trial> MakeTrialReptate(const argtype &args) {
  return MakeTrialMove(std::make_shared<SelectReptate>(args),
    std::make_shared<PerturbReptate>(args),
    "TrialReptate",
    args);
}

std::shared_ptr<Trial> MakeTrialSwapSites(const argtype &args) {
  auto trial = MakeTrial(args);
  trial->set_description("TrialSwapSites");
  trial->set(MakeTrialComputeMove());
  Arguments args_(args);
  args_.dont_check();
  const int site_type1 = args_.key("site_type1").integer();
  const int site_type2 = args_.key("site_type2").integer();
  ASSERT(site_type1 != site_type2, "site types should not match: " <<
    site_type1 << " " << site_type2);
  const std::string part_type = args_.key("particle_type").str();
  trial->add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type1)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type2)}}),
    args
  );
  trial->add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type2)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type1)}}),
    args
  );
  return trial;
}

}  // namespace feasst
