#include <cmath>
#include <vector>
#include "configuration/include/domain.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"

namespace feasst {

Domain::Domain(argtype args) {
  bool is_cubic = used("cubic_box_length", args);
  if (is_cubic) {
    set_cubic(dble("cubic_box_length", &args));
  }

  std::string start("side_length");
  {
    int dim = dimension();
    std::stringstream key;
    key << start << dim;
    while (used(key.str(), args)) {
      ASSERT(!is_cubic, "cubic_box_length argument should not be used in " <<
        "conjunction with side_length arguments");
      add_side_length(dble(key.str(), &args));
      ++dim;
      ASSERT(dim < 1e8, "dim(" << dim << ") is very high. Infinite loop?");
      key.str("");
      key << start << dim;
    }
  }
  set_xy_(dble("xy", &args, 0.0));
  set_xz_(dble("xz", &args, 0.0));
  set_yz_(dble("yz", &args, 0.0));

  for (int dim = 0; dim < dimension(); ++dim) {
    std::stringstream key;
    key << "periodic" << dim;
    if (used(key.str(), args)) {
      if (!boolean(key.str(), &args)) disable(dim);
    }
  }
  check_all_used(args);
}

Domain& Domain::set_cubic(const double box_length) {
  set_side_lengths(Position().set_vector({box_length, box_length, box_length}));
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
}

void Domain::set_side_length(const int dimension, const double length) {
  side_lengths_.set_coord(dimension, length);
}

void Domain::add_side_length(const double length) {
  side_lengths_.push_back(length);
  periodic_.push_back(true);
  resize_opt_(side_lengths_.dimension());
}

void Domain::set_xy_(const double xy) {
  xy_ = xy;
  if (std::abs(xy_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
}

void Domain::set_xz_(const double xz) {
  xz_ = xz;
  if (std::abs(xz_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
}

void Domain::set_yz_(const double yz) {
  yz_ = yz;
  if (std::abs(yz_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
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
  wrap_opt(position, opt_origin_, &opt_rel_, &opt_pbc_, &opt_r2_);
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
  ASSERT(!is_tilted(), "implement triclinic");
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

std::string Domain::status_header() const {
  std::stringstream ss;
  ss << ",volume";
  return ss.str();
}

std::string Domain::status() const {
  std::stringstream ss;
  ss << "," << volume();
  return ss.str();
}

void Domain::serialize(std::ostream& sstr) const {
  feasst_serialize_version(841, sstr);
  feasst_serialize_fstobj(side_lengths_, sstr);
  feasst_serialize(xy_, sstr);
  feasst_serialize(xz_, sstr);
  feasst_serialize(yz_, sstr);
  feasst_serialize(is_tilted_, sstr);
  feasst_serialize(periodic_, sstr);
}

Domain::Domain(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 841, "version mismatch: " << version);
  feasst_deserialize_fstobj(&side_lengths_, sstr);
  feasst_deserialize(&xy_, sstr);
  feasst_deserialize(&xz_, sstr);
  feasst_deserialize(&yz_, sstr);
  feasst_deserialize(&is_tilted_, sstr);
  feasst_deserialize(&periodic_, sstr);
  resize_opt_(side_lengths_.dimension());
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
      (*dbc)[dim] -= side_length*std::rint((*dxv)[dim]/side_length);
      (*dxv)[dim] += (*dbc)[dim];
    }
    *r2 += (*dxv)[dim]*(*dxv)[dim];
  }
}

void Domain::wrap_triclinic_opt(const Position& pos1,
    const Position& pos2,
    Position * rel,
    Position * pbc,
    double * r2) const {
  *r2 = 0;
  const std::vector<double>& side = side_lengths_.coord();
  std::vector<double>* dxv = (*rel).get_coord();
  std::vector<double>* dbc = (*pbc).get_coord();
  if (pos1.dimension() >= 3) {
    (*dxv)[2] = pos1.coord()[2] - pos2.coord()[2];
    if (periodic_[2]) {
      const double side_length = side[2];
      const int num_wrap = std::rint((*dxv)[2]/side_length);
      (*dbc)[2] -= num_wrap*side_length;
      (*dbc)[1] -= num_wrap*yz_;
      (*dbc)[0] -= num_wrap*xz_;
      (*dxv)[2] += (*dbc)[2];
      (*dxv)[1] += (*dbc)[1];
      (*dxv)[0] += (*dbc)[0];
    }
    *r2 += (*dxv)[2]*(*dxv)[2];
  }
  (*dxv)[1] = pos1.coord()[1] - pos2.coord()[1];
  if (periodic_[1]) {
    const double side_length = side[1];
    const int num_wrap = std::rint((*dxv)[1]/side_length);
    (*dbc)[1] -= num_wrap*side_length;
    (*dbc)[0] -= num_wrap*xy_;
    (*dxv)[1] += (*dbc)[1];
    (*dxv)[0] -= (*dbc)[0];
  }
  (*dxv)[0] = pos1.coord()[0] - pos2.coord()[0];
  if (periodic_[0]) {
    const double side_length = side[0];
    (*dbc)[0] -= side_length*std::rint((*dxv)[0]/side_length);
    (*dxv)[0] += (*dbc)[0];
  }
  *r2 += (*dxv)[0]*(*dxv)[0] + (*dxv)[1]*(*dxv)[1];
}

}  // namespace feasst
