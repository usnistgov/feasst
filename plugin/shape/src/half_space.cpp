#include "utils/include/serialize.h"
#include "shape/include/half_space.h"

namespace feasst {

class MapHalfSpace {
 public:
  MapHalfSpace() {
    auto obj = MakeHalfSpace({
      {"dimension", "1"},
      {"intersection", "1"},
      {"direction", "1"},
    });
    obj->deserialize_map()["HalfSpace"] = obj;
  }
};

static MapHalfSpace mapper_ = MapHalfSpace();

HalfSpace::HalfSpace(const argtype &args) : Shape() {
  class_name_ = "HalfSpace";
  args_.init(args);
  dimension_ = args_.key("dimension").integer();
  intersection_ = args_.key("intersection").dble();
  direction_ = args_.key("direction").integer();
  DEBUG("dir " << direction_);
  ASSERT(direction_ == -1 or direction_ == 1, "invalid direction: "
    << direction_);
}

double HalfSpace::nearest_distance(const Position& point) const {
  DEBUG("c " << point.coord(dimension_) << " " <<
                direction_*(intersection_ - point.coord(dimension_)));
  return direction_*(intersection_ - point.coord(dimension_));
}

void HalfSpace::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_half_space_(ostr);
}

void HalfSpace::serialize_half_space_(std::ostream& ostr) const {
  serialize_shape_(ostr);
  feasst_serialize_version(376, ostr);
  feasst_serialize(dimension_, ostr);
  feasst_serialize(intersection_, ostr);
  feasst_serialize(direction_, ostr);
}

HalfSpace::HalfSpace(std::istream& istr) : Shape(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(376 == version, version);
  feasst_deserialize(&dimension_, istr);
  feasst_deserialize(&intersection_, istr);
  feasst_deserialize(&direction_, istr);
}

}  // namespace feasst
