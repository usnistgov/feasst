
#include <algorithm>
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

}  // namespace feasst
