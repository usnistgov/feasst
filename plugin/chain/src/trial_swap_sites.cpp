#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "chain/include/trial_swap_sites.h"

namespace feasst {

FEASST_MAPPER(TrialSwapSites, argtype({{"particle_type", "0"},
  {"site_type1", "0"}, {"site_type2", "1"}}));

TrialSwapSites::TrialSwapSites(argtype * args) : Trial(args) {
  class_name_ = "TrialSwapSites";
  set_description("TrialSwapSites");
  const int site_type1 = integer("site_type1", args);
  const int site_type2 = integer("site_type2", args);
  ASSERT(site_type1 != site_type2, "site types should not match: " <<
    site_type1 << " " << site_type2);
  const std::string part_type = str("particle_type", args);
  argtype stage_args = *args;
  add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type1)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type2)}}),
    &stage_args);
  feasst_check_all_used(stage_args);
  add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type2)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type1)}}),
    args
  );
  set(MakeTrialComputeMove());
}
TrialSwapSites::TrialSwapSites(argtype args) : TrialSwapSites(&args) {
  feasst_check_all_used(args);
}

TrialSwapSites::TrialSwapSites(std::istream& istr) : Trial(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9478, "mismatch version: " << version);
}

void TrialSwapSites::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(9478, ostr);
}

}  // namespace feasst
