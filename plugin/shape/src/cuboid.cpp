#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/cuboid.h"

namespace feasst {

class MapCuboid {
 public:
  MapCuboid() {
    auto obj = MakeCuboid(
      Position({1, 1, 1}),
      Position({0, 0, 0}));
    obj->deserialize_map()["Cuboid"] = obj;
  }
};

static MapCuboid mapper_ = MapCuboid();

Cuboid::Cuboid(const Position& side_lengths,
    const Position& center) : Shape() {
  class_name_ = "Cuboid";
  side_lengths_ = side_lengths;
  center_ = center;
}

double Cuboid::nearest_distance(const Position& point) const {
  FATAL("not implemented");
}

void Cuboid::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_shape_(ostr);
  feasst_serialize_version(6972, ostr);
  feasst_serialize_fstobj(side_lengths_, ostr);
  feasst_serialize_fstobj(center_, ostr);
}

Cuboid::Cuboid(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(6972 == version, version);
  feasst_deserialize_fstobj(&side_lengths_, istr);
  feasst_deserialize_fstobj(&center_, istr);
}

double Cuboid::surface_area() const {
  ASSERT(side_lengths_.dimension() == 3, "assumes 3D");
  return 2*(side_lengths_.coord(0)*side_lengths_.coord(1) +
            side_lengths_.coord(0)*side_lengths_.coord(2) +
            side_lengths_.coord(1)*side_lengths_.coord(2));
}

double Cuboid::volume() const { return product(side_lengths_.coord()); }

}  // namespace feasst
