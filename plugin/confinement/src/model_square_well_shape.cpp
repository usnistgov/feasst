#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/constants.h"
#include "configuration/include/site.h"
#include "configuration/include/model_params.h"
#include "confinement/include/model_square_well_shape.h"

namespace feasst {

class MapModelSquareWellShape {
 public:
  MapModelSquareWellShape() {
    ModelSquareWellShape().deserialize_map()["ModelSquareWellShape"] =
      std::make_shared<ModelSquareWellShape>();
  }
};

static MapModelSquareWellShape map_model_hard_shape_ = MapModelSquareWellShape();

ModelSquareWellShape::ModelSquareWellShape(std::shared_ptr<Shape> shape,
  const argtype& args) : ModelOneBody(), ShapedEntity(shape) {
  class_name_ = "ModelSquareWellShape";
}

void ModelSquareWellShape::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  ShapedEntity::serialize(ostr);
  feasst_serialize_version(6482, ostr);
}

ModelSquareWellShape::ModelSquareWellShape(std::istream& istr)
  : ModelOneBody(istr), ShapedEntity(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6482, "unrecognized verison: " << version);
}

double ModelSquareWellShape::energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) {
  const int type = site.type();
  const double sigma = model_params.select(sigma_index()).value(type);
  const double cutoff = model_params.select(cutoff_index()).value(type);
  const double distance = -shape()->nearest_distance(wrapped_site);
  if (distance <= 0.5*sigma) {
    return NEAR_INFINITY;
  } else if (distance <= cutoff) {
    const double epsilon = model_params.select(epsilon_index()).value(type);
    return -epsilon;
  } else {
    return 0.;
  }
}

}  // namespace feasst
