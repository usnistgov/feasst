#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments_extra.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "shape/include/cuboid.h"

namespace feasst {

FEASST_MAPPER(Cuboid, argtype({{"side_lengths", "0,0,0"}}));

Cuboid::Cuboid(argtype * args) : Shape() {
  class_name_ = "Cuboid";

  DEBUG("parse cubic_side_length");
  if (used("cubic_side_length", *args)) {
    const double len = dble("cubic_side_length", args);
    ASSERT(len > 0, "len: " << len << " must be > 0");
    side_lengths_ = Position({len, len, len});
  } else {
    const std::string side_lengths = str("side_lengths", args, "");
    if (!side_lengths.empty()) {
      side_lengths_ = Position({{"csv", side_lengths}});
    } else {
      WARN("Deprecate Sphere::side_lengths without comma-separated values.");
      side_lengths_ = Position(parse_dimensional(str("side_length", args), args, 4));
    }
  }

  DEBUG("parse center");
  center_.set_to_origin(side_lengths_.size());
  const std::string center = str("center", args, "");
  if (!center.empty()) {
    if (is_found_in(center, ",")) {
      center_ = Position({{"csv", center}});
    } else {
      WARN("Deprecate Sphere::center without comma-separated values.");
      for (int dim = 0; dim < side_lengths_.size(); ++dim) {
        const std::string key = center+str(dim);
        if (used(key, *args)) {
          center_.set_coord(dim, dble(key, args));
        }
      }
    }
  }
}
Cuboid::Cuboid(argtype args) : Cuboid(&args) {
  feasst_check_all_used(args);
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
