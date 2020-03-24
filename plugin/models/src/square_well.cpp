#include "utils/include/serialize.h"
#include "math/include/constants.h"
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

double SquareWell::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const {
  const double& sigma = model_params.mixed_sigma()[type1][type2];
  const double& epsilon = model_params.mixed_epsilon()[type1][type2];
  if (squared_distance <= sigma*sigma) {
    return NEAR_INFINITY;
  }
  return -epsilon;
}

}  // namespace feasst
