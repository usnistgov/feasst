#include "models/include/square_well.h"

namespace feasst {

class MapSquareWell {
 public:
  MapSquareWell() {
    SquareWell().deserialize_map()["SquareWell"] = MakeSquareWell();
  }
};

static MapSquareWell mapper_ = MapSquareWell();

void SquareWell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(553, ostr);
}

SquareWell::SquareWell(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 553, "unrecognized version: " << version);
}

}  // namespace feasst
