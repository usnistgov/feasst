#include <cmath>
#include <vector>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"

namespace feasst {

Domain::Domain(argtype * args) {
  bool is_cubic = used("cubic_side_length", *args);
  if (is_cubic) {
    set_cubic(dble("cubic_side_length", args));
  }
  bool is_cubic_box = used("cubic_box_length", *args);
  ASSERT(!(is_cubic && is_cubic_box),
    "Cannot use both cubic_side_length and cubic_box_length");
  if (is_cubic_box) {
    WARN("cubic_box_length is deprecated. Use cubic_side_length instead.");
    set_cubic(dble("cubic_box_length", args));
  }

  std::string start("side_length");
  {
    if (used(start, *args)) {
      ASSERT(!is_cubic, "cubic_side_length argument should not be used in " <<
        "conjunction with the side_length argument.");
      for (const std::string& st : split(feasst::str(start, args), ',')) {
        add_side_length(str_to_double(st));
      }
    }
    int dim = dimension();
    std::stringstream key;
    key << start << dim;
    while (used(key.str(), *args)) {
      ASSERT(!is_cubic, "cubic_side_length argument should not be used in " <<
        "conjunction with side_length arguments");
      WARN("side_length0 is deprecated. Use comma-separated list for side_length.");
      add_side_length(dble(key.str(), args));
      ++dim;
      ASSERT(dim < 1e8, "dim(" << dim << ") is very high. Infinite loop?");
      key.str("");
      key << start << dim;
    }
  }
  set_xy_(dble("xy", args, 0.0));
  set_xz_(dble("xz", args, 0.0));
  set_yz_(dble("yz", args, 0.0));

  start.assign("periodic");
  if (used(start, *args)) {
    std::vector<std::string> prds = split(feasst::str(start, args), ',');
    ASSERT(dimension() == static_cast<int>(prds.size()), "The number of " <<
      "comma-separated values for periodic argument:" << prds.size() <<
      " must equal the number of dimensions:" << dimension());
    for (int dim = 0; dim < dimension(); ++dim) {
      if (!str_to_bool(prds[dim])) {
        disable(dim);
      }
    }
  } else {
    for (int dim = 0; dim < dimension(); ++dim) {
      std::stringstream key;
      key << start << dim;
      if (used(key.str(), *args)) {
        WARN("periodic[i] is deprecated. Use comma-separated periodic argument"
          << " where the number of values must equal the number of dimensions");
        if (!boolean(key.str(), args)) disable(dim);
      }
    }
  }
  update_h_();
}
Domain::Domain(argtype args) : Domain(&args) {
  feasst_check_all_used(args);
}

Domain& Domain::set_cubic(const double box_length) {
  set_side_lengths(Position().set_vector({box_length, box_length, box_length}));
  update_h_();
  return *this;
}

void Domain::resize_opt_(const int dimension) {
  opt_origin_.set_to_origin(dimension);
  opt_rel_.set_to_origin(dimension);
  opt_pbc_.set_to_origin(dimension);
}

void Domain::set_side_lengths(const Position& side_lengths) {
  side_lengths_ = side_lengths;
  for (int dim = static_cast<int>(periodic_.size());
       dim < static_cast<int>(side_lengths_.size());
       ++dim) {
     periodic_.push_back(true);
  }
  ASSERT(static_cast<int>(periodic_.size()) ==
         static_cast<int>(side_lengths_.size()), "size error");
  resize_opt_(side_lengths_.dimension());
  update_h_();
}

void Domain::set_side_length(const int dimension, const double length) {
  side_lengths_.set_coord(dimension, length);
  update_h_();
}

void Domain::add_side_length(const double length) {
  side_lengths_.push_back(length);
  periodic_.push_back(true);
  resize_opt_(side_lengths_.dimension());
  update_h_();
}

void Domain::set_xy_(const double xy) {
  xy_ = xy;
  if (std::abs(xy_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
  update_h_();
}

void Domain::set_xz_(const double xz) {
  xz_ = xz;
  if (std::abs(xz_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
  update_h_();
}

void Domain::set_yz_(const double yz) {
  yz_ = yz;
  if (std::abs(yz_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
  update_h_();
}

double Domain::volume() const {
  double vol = 1.;
  for (double length : side_lengths_.coord()) {
    vol *= length;
  }
  return vol;
}

Position Domain::shift(const Position& position) const {
  // use the optimized version for consistency
  Position pos2, rel, pbc;
  pos2.set_to_origin_3D();
  pbc.set_to_origin_3D();
  rel = position;
  double r2;
  wrap_opt(position, pos2, &rel, &pbc, &r2);
  rel.subtract(position);
  return rel;
}

const Position& Domain::shift_opt(const Position& position) {
  if (is_tilted_) {
    DEBUG("position " << position.str());
    DEBUG("h_inv " << h_inv_.str());
    h_inv_.multiply(position, &opt_pbc_);
    for (int dim = 0; dim < dimension(); ++dim) {
      const double val = opt_pbc_.coord(dim);
      opt_pbc_.set_coord(dim, val - std::rint(val));
    }
    h_.multiply(opt_pbc_, &opt_rel_);
  } else {
    wrap_opt(position, opt_origin_, &opt_rel_, &opt_pbc_, &opt_r2_);
  }
  opt_rel_.subtract(position);
  return opt_rel_;
}

void Domain::wrap(Position * position) const {
  position->add(shift(*position));
}

Position Domain::random_position(Random * random) const {
  Position position;
  random_position(&position, random);
  return position;
}

void Domain::random_position(Position * position, Random * random) const {
  DEBUG("side_lengths_ " << side_lengths_.str());
  return random->position_in_cuboid(side_lengths_, position);
}

void Domain::unwrap(const int dim, const int num_wrap, Position * shift) const {
  // do nothing if not periodic
  if (!periodic(dim)) {
    return;
  }

  // simple method if cuboid
  if (!is_tilted()) {
    shift->add_to_coord(dim, num_wrap*side_length(dim));
    return;
  }

  // otherwise, unwrap triclinic box
  if (dim == 2) {
    shift->add_to_coord(2, num_wrap*side_length(dim));
    shift->add_to_coord(1, num_wrap*yz());
    shift->add_to_coord(0, num_wrap*xz());
  } else if (dim == 1) {
    shift->add_to_coord(1, num_wrap*side_length(dim));
    shift->add_to_coord(0, num_wrap*xy());
  } else if (dim == 0) {
    shift->add_to_coord(0, num_wrap*side_length(dim));
  } else {
    ERROR("unrecognized dim:" << dim);
  }
}

bool Domain::is_cubic() const {
  if (side_lengths_.size() == 0) {
    return false;
  }
  const double length0 = side_lengths_.coord(0);
  for (const double& len : side_lengths_.coord()) {
    if (std::abs(len - length0) > NEAR_ZERO) {
      return false;
    }
  }
  return true;
}

double Domain::min_side_length() const {
  ASSERT(side_lengths_.dimension() > 0, "no side lengths");
  return minimum(side_lengths_.coord());
}

double Domain::max_side_length() const {
  ASSERT(side_lengths_.dimension() > 0, "no side lengths");
  return maximum(side_lengths_.coord());
}

double Domain::inscribed_sphere_diameter() const {
  double diameter = -1.;
  if (!is_tilted()) {
    diameter = min_side_length();
  } else if (dimension() == 3) {
    Position A({side_length(0), 0., 0.}),
             B({xy(), side_length(1), 0.}),
             C({xz(), yz(), side_length(2)});
    DEBUG("A " << A.str());
    DEBUG("B " << B.str());
    DEBUG("C " << C.str());
    Position BcrossC = B.cross_product(C);
    Position CcrossA = C.cross_product(A);
    Position AcrossB = A.cross_product(B);
    const double widthA = BcrossC.dot_product(A)/BcrossC.distance();
    DEBUG("widthA " << widthA);
    const double widthB = CcrossA.dot_product(B)/CcrossA.distance();
    DEBUG("widthB " << widthB);
    const double widthC = AcrossB.dot_product(C)/AcrossB.distance();
    DEBUG("widthC " << widthC);
    if (widthA <= widthB && widthA <= widthC) {
      diameter = widthA;
    } else if (widthB <= widthA && widthB <= widthC) {
      diameter = widthB;
    } else {
      diameter = widthC;
    }
  } else if (dimension() == 2) {
    const double ly = side_length(1);
    diameter = side_length(0)*ly/std::sqrt(xy()*xy() + ly*ly);
    if (diameter > ly) {
      diameter = ly;
    }
  } else {
    FATAL("dimension: " << dimension() << " not implemented");
  }
  return diameter;
}

std::string Domain::status_header(const std::string append) const {
  std::stringstream ss;
  ss << ",volume" << append;
  return ss.str();
}

std::string Domain::status() const {
  std::stringstream ss;
  ss << "," << volume();
  return ss.str();
}

void Domain::serialize(std::ostream& sstr) const {
  feasst_serialize_version(842, sstr);
  feasst_serialize_fstobj(side_lengths_, sstr);
  feasst_serialize(xy_, sstr);
  feasst_serialize(xz_, sstr);
  feasst_serialize(yz_, sstr);
  feasst_serialize(is_tilted_, sstr);
  feasst_serialize(periodic_, sstr);
}

Domain::Domain(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version >= 841 && version <= 842, "version mismatch: " << version);
  feasst_deserialize_fstobj(&side_lengths_, sstr);
  feasst_deserialize(&xy_, sstr);
  feasst_deserialize(&xz_, sstr);
  feasst_deserialize(&yz_, sstr);
  feasst_deserialize(&is_tilted_, sstr);
  feasst_deserialize(&periodic_, sstr);
  resize_opt_(side_lengths_.dimension());
  update_h_();
}

void Domain::wrap_opt(const Position& pos1,
    const Position& pos2,
    Position * rel,
    Position * pbc,
    double * r2) const {
  if (is_tilted_) {
    wrap_triclinic_opt(pos1, pos2, rel, pbc, r2);
    return;
  }
  const int dimen = pos1.dimension();
  *r2 = 0;
  const std::vector<double>& side = side_lengths_.coord();
  std::vector<double>* dxv = (*rel).get_coord();
  std::vector<double>* dbc = (*pbc).get_coord();
  for (int dim = 0; dim < dimen; ++dim) {
    (*dxv)[dim] = pos1.coord()[dim] - pos2.coord()[dim];
    (*dbc)[dim] = 0.;
    const double side_length = side[dim];
    if (periodic_[dim]) {
      const double dx = side_length*std::rint((*dxv)[dim]/side_length);
      (*dbc)[dim] -= dx;
      (*dxv)[dim] -= dx;
    }
    const double dxvdim = (*dxv)[dim];
    *r2 += dxvdim*dxvdim;
  }
}

void Domain::wrap_triclinic_opt(const Position& pos1,
    const Position& pos2,
    Position * rel,
    Position * pbc,
    double * r2) const {
  DEBUG("wrapping triclinc opt " << pos1.str() << " " << pos2.str());
  *r2 = 0;
  const std::vector<double>& side = side_lengths_.coord();
  std::vector<double>* dxv = (*rel).get_coord();
  std::vector<double>* dbc = (*pbc).get_coord();
  (*dbc)[0] = 0.;
  (*dbc)[1] = 0.;
  (*dxv)[0] = pos1.coord()[0] - pos2.coord()[0];
  (*dxv)[1] = pos1.coord()[1] - pos2.coord()[1];
  if (pos1.dimension() >= 3) {
    (*dxv)[2] = pos1.coord()[2] - pos2.coord()[2];
    if (periodic_[2]) {
      (*dbc)[2] = 0.;
      const double side_length = side[2];
      const int num_wrap = std::rint((*dxv)[2]/side_length);
      const double dz = num_wrap*side_length;
      const double dy = num_wrap*yz_;
      const double dx = num_wrap*xz_;
      (*dbc)[2] -= dz;
      (*dbc)[1] -= dy;
      (*dbc)[0] -= dx;
      (*dxv)[2] -= dz;
      (*dxv)[1] -= dy;
      (*dxv)[0] -= dx;
    }
    const double dxv2 = (*dxv)[2];
    *r2 += dxv2*dxv2;
  }
  if (periodic_[1]) {
    const double side_length = side[1];
    const int num_wrap = std::rint((*dxv)[1]/side_length);
    const double dy = num_wrap*side_length;
    const double dx = num_wrap*xy_;
    (*dbc)[1] -= dy;
    (*dbc)[0] -= dx;
    (*dxv)[1] -= dy;
    (*dxv)[0] -= dx;
  }
  if (periodic_[0]) {
    const double side_length = side[0];
    const double dx = side_length*std::rint((*dxv)[0]/side_length);
    (*dbc)[0] -= dx;
    (*dxv)[0] -= dx;
  }
  const double dxv0 = (*dxv)[0];
  const double dxv1 = (*dxv)[1];
  *r2 += dxv0*dxv0 + dxv1*dxv1;
}

void Domain::update_h_() {
  DEBUG("dim " << dimension());
  if (dimension() == 2) {
    const double lx = side_lengths_.coord(0);
    const double ly = side_lengths_.coord(1);
    h_.set_size(2, 2);
    h_.set_value(0, 0, lx);
    h_.set_value(0, 1, xy_);
    h_.set_value(1, 0, 0.);
    h_.set_value(1, 1, ly);
    h_inv_.set_size(2, 2);
    h_inv_.set_value(0, 0, 1./lx);
    h_inv_.set_value(0, 1, -xy_/lx/ly);
    h_inv_.set_value(1, 0, 0.);
    h_inv_.set_value(1, 1, 1./ly);
  } else if (dimension() == 3) {
    const double lx = side_lengths_.coord(0);
    const double ly = side_lengths_.coord(1);
    const double lz = side_lengths_.coord(2);
    h_.set_size(3, 3);
    h_.set_value(0, 0, lx);
    h_.set_value(0, 1, xy_);
    h_.set_value(0, 2, xz_);
    h_.set_value(1, 0, 0.);
    h_.set_value(1, 1, ly);
    h_.set_value(1, 2, yz_);
    h_.set_value(2, 0, 0.);
    h_.set_value(2, 1, 0.);
    h_.set_value(2, 2, lz);
    h_inv_.set_size(3, 3);
    h_inv_.set_value(0, 0, 1./lx);
    h_inv_.set_value(0, 1, -xy_/lx/ly);
    h_inv_.set_value(0, 2, (xy_*yz_ - xz_*ly)/lx/ly/lz);
    h_inv_.set_value(1, 0, 0.);
    h_inv_.set_value(1, 1, 1./ly);
    h_inv_.set_value(1, 2, -yz_/ly/lz);
    h_inv_.set_value(2, 0, 0.);
    h_inv_.set_value(2, 1, 0.);
    h_inv_.set_value(2, 2, 1./lz);
    DEBUG("h rows " << h_.num_rows());
    DEBUG("h columns " << h_.num_columns());
  }
  DEBUG("h " << h_.str());
  DEBUG("h_inv " << h_inv_.str());
}

void Domain::cartesian2scaled(const Position& cartesian, Position * scaled) const {
  h_inv_.multiply(cartesian, scaled);
}

void Domain::cartesian2scaled_wrap(const Position& cartesian, Position * scaled) const {
  cartesian2scaled(cartesian, scaled);
  for (int dim = 0; dim < dimension(); ++dim) {
    const double val = scaled->coord(dim);
    scaled->set_coord(dim, val - std::rint(val));
  }
}

bool Domain::is_minimum_image_for_cutoff(const double cutoff) const {
  DEBUG("checking if cutoff " << cutoff << " follows minimum image.");
  if (is_tilted()) {
    // titled domains do not check for periodicity
    if (cutoff - NEAR_ZERO > 0.5*inscribed_sphere_diameter()) {
      return false;
    }
  } else {
    for (int dim = 0; dim < dimension(); ++dim) {
      DEBUG("dim " << dim << " periodic " << periodic(dim));
      if (periodic(dim)) {
        if (cutoff - NEAR_ZERO > 0.5*side_length(dim)) {
          return false;
        }
      }
    }
  }
  return true;
}

void Domain::wrap_opt(const Position& unwrapped, Position * wrapped, Position * scaled) const {
  for (Position * pos : {wrapped, scaled}) {
    if (pos->size() != unwrapped.size()) {
      pos->set_to_origin(unwrapped.size());
    }
  }
  cartesian2scaled_wrap(unwrapped, scaled);
  h_.multiply(*scaled, wrapped);
}

}  // namespace feasst
