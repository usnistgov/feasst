#include "system/include/lennard_jones.h"

namespace feasst {

class MapLennardJones {
 public:
  MapLennardJones() {
    LennardJones().deserialize_map()["LennardJones"] = MakeLennardJones();
  }
};

static MapLennardJones mapper_ = MapLennardJones();

void LennardJones::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_lennard_jones_(ostr);
}

void LennardJones::serialize_lennard_jones_(std::ostream& ostr) const {
  feasst_serialize_version(763, ostr);
  feasst_serialize(hard_sphere_threshold_sq_, ostr);
}

LennardJones::LennardJones(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(763 == version, version);
  feasst_deserialize(&hard_sphere_threshold_sq_, istr);
}

}  // namespace feasst
