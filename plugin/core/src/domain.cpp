
#include "core/include/domain.h"
#include "core/include/constants.h"

namespace feasst {

double Domain::volume() const {
  double vol = 1.;
  //for (int dim = 0; dim < static_cast<int>(side_length_.size()); ++dim) {
    //vol *= side_length_.coord(dim);
  for (double length : side_length_.coord()) {
    vol *= length;
  }
  return vol;
}

void Domain::wrap(Position * position) const {
  position->add(shift(*position));
}

double Domain::min_side_length() const {
  double min = NEAR_INFINITY;
  for (double coord : side_length_.coord()) {
    if (coord < min) {
      min = coord;
    }
  }
  return min;
}

int Domain::cell_id(const Position& position,
                    const Cells& cells) const {
  Position scaled(position);
  wrap(&scaled);
  scaled.divide(side_length());
  return cells.id(scaled.coord());
}

void Domain::init_cells(const double min_length,
                        const int group_index) {
  ASSERT(group_index >= 0, "error");
  Cells cell;
  cell.create(min_length, side_length().coord());
  std::stringstream ss;
  ss << "cell";
  ss << cells_.size();
  cell.set_label(ss.str());
  cell.add_property("group", group_index);
  cells_.push_back(cell);
}

}  // namespace feasst
