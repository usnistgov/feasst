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
    const argtype& args) {
  Arguments args_(args);

  // read alpha and epsilon
  std::vector<double> alpha, epsilon;
  std::string start;
  start.assign("alpha");
  if (args_.key(start).used()) {
    alpha.push_back(args_.dble());
    epsilon.push_back(args_.key("epsilon").dble());
  } else {
    int type = static_cast<int>(alpha.size());
    std::stringstream key;
    key << start << type;
    while (args_.key(key.str()).used()) {
      alpha.push_back(args_.dble());
      key.str("");
      key << "epsilon" << type;
      epsilon.push_back(args_.key(key.str()).dble());
      ++type;
      key.str("");
      key << start << type;
    }
  }

  const bool invert = args_.key("invert").dflt("true").boolean();
  const double max_radius = args_.key("max_radius").dble();
  const int num_radius = args_.key("num_radius").integer();
  const double density = args_.key("density").dble();
  double sum = 0.;
  ASSERT(point.dimension() == 3, "assumes 3d");
  const double dr = max_radius/static_cast<double>(num_radius);
  Position surf_point(3);
  double lower_radius = 0.;
  double lower_vol = 0.;
  double lower_frac = 0.;
  bool ins = is_inside(point);
  //if ((!invert && ins) || (invert && !ins)) lower_frac = 1.;
  if (!spheres_ || !meshes_) {
    spheres_ = std::make_shared<std::vector<std::shared_ptr<Sphere> > >();
    meshes_ = std::make_shared<std::vector<std::vector<Position> > >();
    for (double radius = dr; radius < max_radius + dr/2.; radius += dr) {
      spheres_->push_back(MakeSphere({{"radius", str(radius)}}));
      meshes_->push_back(std::vector<Position>());
      spheres_->back()->surface_mesh(density, &meshes_->back());
    }
  }
  int irad = 0;
  for (double radius = dr; radius < max_radius + dr/2.; radius += dr) {
    int inside = 0;
    DEBUG("irad " << irad);
    DEBUG("radius " << radius);
    DEBUG("lower_frac " << lower_frac);
    DEBUG("lower_radius " << lower_radius);
    DEBUG("lower_vol " << lower_vol);
    int num_points = static_cast<int>((*meshes_)[irad].size());
    DEBUG("num_points " << num_points);
    if (num_points > 0) {
      const RotationMatrix rot_mat = random->rotation(3, 180);
      for (int ip = 0; ip < num_points; ++ip) {
        rot_mat.multiply((*meshes_)[irad][ip], &surf_point);
        surf_point.add(point);
        TRACE(surf_point.str());
        ins = is_inside(surf_point);
        if ((!invert && ins) || (invert && !ins)) ++inside;
      }
      const double frac = static_cast<double>(inside)/static_cast<double>(num_points);
      const double volume = (*spheres_)[irad]->volume();
      DEBUG("frac " << frac);
      DEBUG("vol " << volume);
      const double shell_vol = volume - lower_vol;
      DEBUG("shell_vol " << shell_vol);
      const double radius_av = 0.5*(radius + lower_radius);
      double frac_av = 0.5*(frac + lower_frac);
      if (irad == 0) frac_av = frac;
      DEBUG("frac_av " << frac_av);
      for (int ith = 0; ith < static_cast<int>(alpha.size()); ++ith) {
        sum += epsilon[ith]*frac_av*shell_vol*std::pow(radius_av, -alpha[ith]);
      }
      DEBUG("sum: " << sum);
      lower_radius = radius;
      lower_vol = volume;
      lower_frac = frac;
    }
    ++irad;
  }
  return sum;
}

}  // namespace feasst
