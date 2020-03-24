#include <algorithm>
#include "utils/include/serialize.h"
#include "growth_expanded/include/macrostate_growth_expanded.h"

namespace feasst {

class MapMacrostateGrowthExpanded {
 public:
  MapMacrostateGrowthExpanded() {
    MacrostateGrowthExpanded(Histogram()).deserialize_map()["MacrostateGrowthExpanded"] =
      MakeMacrostateGrowthExpanded(Histogram());
  }
};

static MapMacrostateGrowthExpanded mapper_ = MapMacrostateGrowthExpanded();

std::shared_ptr<Macrostate> MacrostateGrowthExpanded::create(std::istream& istr) const {
  return std::make_shared<MacrostateGrowthExpanded>(istr);
}

MacrostateGrowthExpanded::MacrostateGrowthExpanded(std::istream& istr)
  : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6030, "unrecognized version: " << version);
}

void MacrostateGrowthExpanded::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(6030, ostr);
}

double MacrostateGrowthExpanded::value(const System* system,
    const Criteria* criteria) {
  const double num_particles = static_cast<double>(system->configuration().num_particles());
  double rtrn;
  DEBUG("nump " << num_particles);
  DEBUG("state " << criteria->trial_state());
  DEBUG("num state: " << criteria->num_trial_states());
  if (criteria->trial_state() == 0) {
    rtrn = num_particles;
  } else {
    rtrn = num_particles - 1 +
        (static_cast<double>(criteria->trial_state())/
         static_cast<double>(criteria->num_trial_states()));
  }
  DEBUG("rtrn: " << rtrn);
  return rtrn;
}

}  // namespace feasst
