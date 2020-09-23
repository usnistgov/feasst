#include <cmath>
#include <vector>
#include "configuration/include/domain.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"

namespace feasst {

Domain::Domain(const argtype& args) {
  Arguments args_(args);

  bool is_cubic = args_.key("cubic_box_length").used();
  if (is_cubic) {
    set_cubic(args_.dble());
  }

  std::string start("side_length");
  {
    int dim = dimension();
    std::stringstream key;
    key << start << dim;
    while (args_.key(key.str()).used()) {
      ASSERT(!is_cubic, "cubic_box_length argument should not be used in " <<
        "conjunction with side_length arguments");
      add_side_length(args_.dble());
      ++dim;
      ASSERT(dim < 1e8, "dim(" << dim << ") is very high. Infinite loop?");
      key.str("");
      key << start << dim;
    }
  }
  set_xy_(args_.key("xy").dflt("0.0").dble());
  set_xz_(args_.key("xz").dflt("0.0").dble());
  set_yz_(args_.key("yz").dflt("0.0").dble());

  for (int dim = 0; dim < dimension(); ++dim) {
    std::stringstream key;
    key << "periodic" << dim;
    if (args_.key(key.str()).used()) {
      if (!args_.boolean()) disable(dim);
    }
  }

  DEBUG("parse cells");
  // HWH make this more modular
  start.assign("init_cells");
  // if only one cell, drop subscript
  if (args_.key(start).used()) {
    const double min_length = args_.dble();
    int group_index = args_.key("cell_group").dflt("0").integer();
    init_cells(min_length, group_index);
  } else {
    int cell = 0;
    std::stringstream key;
    key << start << cell;
    while (args_.key(key.str()).used()) {
      const double min_length = args_.dble();
      std::stringstream cgrp;
      cgrp << "cell_group" << cell;
      int group_index = args_.key(cgrp.str()).dflt("0").integer();
      init_cells(min_length, group_index);
      ++cell;
      ASSERT(cell < 1e8, "cell(" << cell << ") is very high. Infinite loop?");
      key.str("");
      key << start << cell;
    }
  }
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
  ASSERT(cells_.size() == 0, "Implement check of cell list when changing box");
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

// HWH note if there are problems with scaled coordinates here, it probably
// means there is an issue with wrapping. As currently implemented, translations
// automatically wrap. So if you're doing a test without them you might run
// into this issue.
int Domain::cell_id(const Position& position,
                    const Cells& cells) {
//  Position scaled(position);
//  DEBUG("scaled before wrap " << scaled.str() << " pos " << position.str() <<
//    " box " << side_lengths().str());
  wrap_opt(position, opt_origin_, &opt_rel_, &opt_pbc_, &opt_r2_);
  //wrap(&scaled);
  DEBUG("opt_rel_ after wrap " << opt_rel_.str() << " pos " << position.str() <<
    " box " << side_lengths().str());
  opt_rel_.divide(side_lengths());
  DEBUG("opt_rel_ " << opt_rel_.str() << " pos " << position.str() << " box "
    << side_lengths().str());
  return cells.id(opt_rel_.coord());
}

void Domain::init_cells(const double min_length,
                        const int group_index) {
  ASSERT(side_lengths().size() > 0, "cannot define cells before domain sides");
  ASSERT(!is_tilted(), "implement triclinic");
  ASSERT(group_index >= 0, "error");
  Cells cell;
  cell.create(min_length, side_lengths().coord());
  cell.set_type(cells_.size());
  cell.set_group(group_index);
  if (cell.num_total() > 0) {
    cells_.push_back(cell);
  } else {
    INFO("Requested cell list rejected: min_length:" << min_length <<
         " did not meet requirements.");
  }
}

const Cells& Domain::cells(const int index) const {
  ASSERT(index < static_cast<int>(cells_.size()), "index error");
  return cells_[index];
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
  feasst_serialize_fstobj(cells_, sstr);
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
  feasst_deserialize_fstobj(&cells_, sstr);
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
