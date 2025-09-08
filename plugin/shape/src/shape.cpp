#include <cmath>
#include "utils/include/serialize_extra.h"
#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
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

void Shape::serialize_shape_(std::ostream& ostr) const {
  feasst_serialize_version(6874, ostr);
}

Shape::Shape(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(6874 == version, "mismatch version: " << version);
}

std::shared_ptr<Shape> Shape::create(std::istream& istr) const {
  FATAL("not implemented"); }
std::shared_ptr<Shape> Shape::create(argtype * args) const {
  FATAL("not implemented"); }

std::shared_ptr<Shape> Shape::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
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
    argtype * args) {

  // read alpha and epsilon
  std::vector<double> alpha, epsilon;
  std::string start;
  start.assign("wall_alpha");
  if (used(start, *args)) {
    alpha.push_back(dble(start, args));
    epsilon.push_back(dble("wall_epsilon", args));
  } else {
    int type = static_cast<int>(alpha.size());
    std::stringstream key;
    key << start << type;
    while (used(key.str(), *args)) {
      alpha.push_back(dble(key.str(), args));
      key.str("");
      key << "wall_epsilon" << type;
      epsilon.push_back(dble(key.str(), args));
      ++type;
      key.str("");
      key << start << type;
    }
  }

  const bool invert = boolean("invert", args, true);
  const double max_radius = dble("max_radius", args);
  const int num_shells = integer("num_shells", args);
  const double points_per_shell = dble("points_per_shell", args);
  double sum = 0.;
  ASSERT(point.dimension() == 3, "assumes 3d");
  const double dr = max_radius/static_cast<double>(num_shells);
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
      spheres_->back()->surface_mesh(points_per_shell, &meshes_->back());
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
  //feasst_check_all_used(*args);
  return sum;
}
double Shape::integrate(
    const Position& point,
    Random * random,
    argtype args) {
  return integrate(point, random, &args);
}

std::vector<Position> Shape::grid(const Position& upper, const Position& lower,
    const int num) {
  std::vector<Position> grid;
  ASSERT(upper.dimension() == 3, "implemented for 3D");
  Position x(3);
  std::vector<double> dx(3);
  for (int dim = 0; dim < 3; ++dim) {
    dx[dim] = (upper.coord(dim) - lower.coord(dim))/(num - 1);
  }
  for ((*x.get_coord())[0] = lower.coord(0);
       x.coord(0) <= upper.coord(0);
       (*x.get_coord())[0] += dx[0]) {
    for ((*x.get_coord())[1] = lower.coord(1);
         x.coord(1) <= upper.coord(1);
         (*x.get_coord())[1] += dx[1]) {
      for ((*x.get_coord())[2] = lower.coord(2);
           x.coord(2) <= upper.coord(2);
           (*x.get_coord())[2] += dx[2]) {
        if (is_inside(x)) {
          grid.push_back(x);
        }
      }
    }
  }
  return grid;
}

std::shared_ptr<Shape> Shape::factory(const std::string name, argtype * args) {
  DEBUG("name: " << name << ", args: " << str(*args));
  return template_factory(deserialize_map(), name, args);
}

}  // namespace feasst
