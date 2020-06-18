#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/random.h"
#include "math/include/matrix.h"
#include "shape/include/shape.h"
#include "shape/include/sphere.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Shape> >& Shape::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Shape> >* ans =
     new std::map<std::string, std::shared_ptr<Shape> >();
  return *ans;
}

bool Shape::is_inside(const Position& point) const {
  TRACE(point.str());
  TRACE(nearest_distance(point));
  if (nearest_distance(point) < 0) {
    return true;
  }
  return false;
}

bool Shape::is_inside(const Position& point, const double diameter) const {
  if (nearest_distance(point) + 0.5*diameter < 0) {
    return true;
  }
  return false;
}

void Shape::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Shape> Shape::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Shape> Shape::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr);
}

void ShapedEntity::serialize(std::ostream& ostr) const {
  feasst_serialize_version(9249, ostr);
  feasst_serialize_fstdr(shape_, ostr);
}

ShapedEntity::ShapedEntity(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9249, "unrecognized verison: " << version);
  // feasst_deserialize_fstdr(shape_, istr);
  { // HWH for unknown reasons the above template function does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      shape_ = shape_->deserialize(istr);
    }
  }
}

double Shape::surface_area() const { FATAL("not implemented"); }

double Shape::volume() const { FATAL("not implemented"); }

double Shape::integrate(
    const Position& point,
    Random * random,
    const argtype& args) const {
  Arguments args_(args);
  const bool invert = args_.key("invert").dflt("true").boolean();
  const double alpha = args_.key("alpha").dflt("6").dble();
  const double max_radius = args_.key("max_radius").dble();
  const int num_radius = args_.key("num_radius").integer();
  const double density = args_.key("density").dble();
  double sum = 0.;
  ASSERT(point.dimension() == 3, "assumes 3d");
  const double dr = max_radius/static_cast<double>(num_radius);
  Position surf_point;
  double lower_radius = 0.;
  double lower_vol = 0.;
  double lower_frac = 0.;
  bool ins = is_inside(point);
  if ((!invert && ins) || (invert && !ins)) lower_frac = 1.;
  std::vector<Position> points;
  for (double radius = dr; radius < max_radius + dr/2.; radius += dr) {
    int inside = 0;
    TRACE("radius " << radius);
    auto sphere = MakeSphere({{"radius", str(radius)}});
    sphere->surface_mesh(density, &points);
    int num_points = static_cast<int>(points.size());
    const RotationMatrix rot_mat = random->rotation(3, 180);
    for (int ip = 0; ip < num_points; ++ip) {
      surf_point = rot_mat.multiply(points[ip]);
      surf_point.add(point);
      TRACE(surf_point.str());
      ins = is_inside(surf_point);
      if ((!invert && ins) || (invert && !ins)) ++inside;
    }
    const double frac = static_cast<double>(inside)/static_cast<double>(num_points);
    const double volume = sphere->volume();
    const double shell_vol = volume - lower_vol;
    const double radius_av = 0.5*(radius + lower_radius);
    sum += 0.5*(frac + lower_frac)*shell_vol*std::pow(radius_av, -alpha);
    lower_radius = radius;
    lower_vol = volume;
    lower_frac = frac;
  }
  return sum;
}


}  // namespace feasst
