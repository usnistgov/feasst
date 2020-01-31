
#include <vector>
#include "configuration/include/domain.h"
#include "math/include/constants.h"

namespace feasst {

Domain::Domain() {
  set_xy();
  set_xz();
  set_yz();
}

Domain& Domain::set_cubic(const double box_length) {
  std::vector<double> cube = {box_length, box_length, box_length};
  Position side_length;
  side_length.set_vector(cube);
  set_side_length(side_length);
  return *this;
}

void Domain::set_side_length(const Position& side_length) {
  side_length_ = side_length;
  for (int dim = static_cast<int>(periodic_.size());
       dim < static_cast<int>(side_length_.size());
       ++dim) {
     periodic_.push_back(true);
  }
  ASSERT(static_cast<int>(periodic_.size()) ==
         static_cast<int>(side_length_.size()), "size error");
}

Domain& Domain::set_xy(const double xy) {
  xy_ = xy;
  if (std::abs(xy_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
  return *this;
}

Domain& Domain::set_xz(const double xz) {
  xz_ = xz;
  if (std::abs(xz_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
  return *this;
}

Domain& Domain::set_yz(const double yz) {
  yz_ = yz;
  if (std::abs(yz_) > NEAR_ZERO) {
    is_tilted_ = true;
  }
  return *this;
}

double Domain::volume() const {
  double vol = 1.;
  for (double length : side_length_.coord()) {
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

void Domain::wrap(Position * position) const {
  position->add(shift(*position));
}

Position Domain::random_position(Random * random) const {
  Position position;
  random_position(&position, random);
  return position;
}

void Domain::random_position(Position * position, Random * random) const {
  DEBUG("side_length_ " << side_length_.str());
  ASSERT(!is_tilted(), "implement triclinic");
  return random->position_in_cuboid(side_length_, position);
}

// HWH note if there are problems with scaled coordinates here, it probably
// means there is an issue with wrapping. As currently implemented, translations
// automatically wrap. So if you're doing a test without them you might run
// into this issue.
int Domain::cell_id(const Position& position,
                    const Cells& cells) const {
  Position scaled(position);
  DEBUG("scaled before wrap " << scaled.str() << " pos " << position.str() << " box " << side_length().str());
  wrap(&scaled);
  DEBUG("scaled after wrap " << scaled.str() << " pos " << position.str() << " box " << side_length().str());
  scaled.divide(side_length());
  DEBUG("scaled " << scaled.str() << " pos " << position.str() << " box " << side_length().str());
  return cells.id(scaled.coord());
}

void Domain::init_cells(const double min_length,
                        const int group_index) {
  ASSERT(side_length().size() > 0, "cannot define cells before domain sides");
  ASSERT(!is_tilted(), "implement triclinic");
  ASSERT(group_index >= 0, "error");
  Cells cell;
  cell.create(min_length, side_length().coord());
  std::stringstream ss;
  ss << "cell";
  ss << cells_.size();
  cell.set_label(ss.str());
  cell.add_property("group", group_index);
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

void Domain::serialize(std::ostream& sstr) const {
  feasst_serialize_version(1, sstr);
  side_length_.serialize(sstr);
  feasst_serialize(xy_, sstr);
  feasst_serialize(xz_, sstr);
  feasst_serialize(yz_, sstr);
  feasst_serialize(is_tilted_, sstr);
  feasst_serialize(periodic_, sstr);
  feasst_serialize_fstobj(cells_, sstr);
}

Domain::Domain(std::istream& sstr) {
  feasst_deserialize_version(sstr);
  side_length_ = Position(sstr);
  feasst_deserialize(&xy_, sstr);
  feasst_deserialize(&xz_, sstr);
  feasst_deserialize(&yz_, sstr);
  feasst_deserialize(&is_tilted_, sstr);
  feasst_deserialize(&periodic_, sstr);
  feasst_deserialize_fstobj(&cells_, sstr);
}

}  // namespace feasst
