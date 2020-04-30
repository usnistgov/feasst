#include "chain/include/trial_swap_sites.h"
#include "utils/include/serialize.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

class MapTrialSwapSites {
 public:
  MapTrialSwapSites() {
    auto obj = MakeTrialSwapSites({{"particle_type", "0"}, {"site_type1", "0"}, {"site_type2", "1"}});
    obj->deserialize_map()["TrialSwapSites"] = obj;
  }
};

static MapTrialSwapSites mapper_ = MapTrialSwapSites();

std::shared_ptr<Trial> TrialSwapSites::create(std::istream& istr) const {
  return std::make_shared<TrialSwapSites>(istr);
}

TrialSwapSites::TrialSwapSites(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialSwapSites", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(346 == version, "mismatch version: " << version);
}

void TrialSwapSites::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(346, ostr);
}

TrialSwapSites::TrialSwapSites(const argtype& args) : Trial(args) {
  class_name_ = "TrialSwapSites";
  set(MakeTrialComputeMove());
  Arguments args_(args);
  args_.dont_check();
  const int site_type1 = args_.key("site_type1").integer();
  const int site_type2 = args_.key("site_type2").integer();
  ASSERT(site_type1 != site_type2, "site types should not match: " <<
    site_type1 << " " << site_type2);
  const std::string part_type = args_.key("particle_type").str();
  add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type1)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type2)}}),
    args
  );
  add_stage(
    MakeSelectSiteOfType({{"site_type", str(site_type2)}, {"particle_type", part_type}}),
    MakePerturbSiteType({{"type", str(site_type1)}}),
    args
  );
}

}  // namespace feasst
