
#include <cmath>
#include <sstream>
#include "configuration/include/cells.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"
#include "utils/include/utils_io.h"

namespace feasst {

void Cells::create(const double min_length,
                   const std::vector<double> side_lengths) {
  ASSERT(min_length > 1e-15, "min_length(" << min_length << ") too small");
  clear();
  for (double side_length : side_lengths) {
    num_.push_back(static_cast<int>(side_length/min_length));
  }
  if (num_total() <= std::pow(3, side_lengths.size())) {
    clear();
    return;
  }
  ASSERT(num_total() < 1e8, "too many cells");
  neighbor_.clear();
  neighbor_.resize(num_total());
  if (static_cast<int>(side_lengths.size()) == 3) {
    build_neighbors_3D_();
  } else if (static_cast<int>(side_lengths.size()) == 2) {
    build_neighbors_2D_();
  } else {
    ASSERT(false, "unrecognized dimension(" << side_lengths.size() << ")");
  }
  build_particles_();
}

void Cells::build_neighbors_3D_() {
  for (int xcell1 = 0; xcell1 < num(0); ++xcell1) {
  for (int ycell1 = 0; ycell1 < num(1); ++ycell1) {
  for (int zcell1 = 0; zcell1 < num(2); ++zcell1) {
    const int cell = id_({xcell1, ycell1, zcell1});
    for (int xcell2 = xcell1 - 1; xcell2 <= xcell1 + 1; ++xcell2) {
    for (int ycell2 = ycell1 - 1; ycell2 <= ycell1 + 1; ++ycell2) {
    for (int zcell2 = zcell1 - 1; zcell2 <= zcell1 + 1; ++zcell2) {
      neighbor_[cell].push_back(id_({xcell2, ycell2, zcell2}));
    }}}
  }}}
}

void Cells::build_neighbors_2D_() {
  for (int xcell1 = 0; xcell1 < num(0); ++xcell1) {
  for (int ycell1 = 0; ycell1 < num(1); ++ycell1) {
    const int cell = id_({xcell1, ycell1});
    for (int xcell2 = xcell1 - 1; xcell2 <= xcell1 + 1; ++xcell2) {
    for (int ycell2 = ycell1 - 1; ycell2 <= ycell1 + 1; ++ycell2) {
      neighbor_[cell].push_back(id_({xcell2, ycell2}));
    }}
  }}
}

void Cells::build_particles_() {
  particles_.resize(num_total());
}

int Cells::num_total() const {
  return product(num_);
}

void Cells::clear() {
  num_.clear();
  neighbor_.clear();
}

int Cells::id_(std::vector<int> position) {
  const int mx = num_[0];
  if (static_cast<int>(position.size()) == 2) {
    const int my = num_[1];
    return (position[0] + mx)%mx +
           mx*((position[1] + my)%my);
  } else if (static_cast<int>(position.size()) == 3) {
    const int my = num_[1];
    const int mz = num_[2];
    return (position[0] + mx )%mx +
           mx*((position[1] + my)%my +
           my*((position[2] + mz)%mz));
  }
  ASSERT(0, "unrecognized dimensionality");
  return -1;
}

int Cells::id(const std::vector<double>& scaled_coord) const {
  ASSERT(scaled_coord.size() == num_.size(), "size error");
  int cell = 0;
  std::vector<int> cells(scaled_coord.size());
  for (int dim = 0; dim < static_cast<int>(cells.size()); ++dim) {
    ASSERT(std::abs(scaled_coord[dim]) <= 0.5,
      MAX_PRECISION << scaled_coord[dim] << " is not "
      << "scaled coordinates");
    cells[dim] = static_cast<int>(
      num_[dim]*(scaled_coord[dim] + 0.5)) % num_[dim];
    double prod = cells[dim];
    for (int dim2 = 0; dim2 < dim; ++dim2) {
      prod *= num_[dim2];
    }
    cell += prod;
  }
  ASSERT( (cell >= 0) and (cell < num_total()),
    "cell(" << cell << ")");
  return cell;
}

int Cells::num_sites() const {
  int num = 0;
  for (const Select& sel : particles()) {
    num += sel.num_sites();
  }
  return num;
}

std::string Cells::str() const {
  std::stringstream ss;
  for (int cell = 0; cell < num_total(); ++cell) {
    const Select& select = particles_[cell];
    if (select.num_particles() != 0) {
      ss << "cell: " << cell << " p: " << select.str();
    }
  }
  return ss.str();
}

void Cells::serialize(std::ostream& sstr) const {
  feasst_serialize_version(1, sstr);
  feasst_serialize(label_, sstr);
  feasst_serialize(num_, sstr);
  feasst_serialize(neighbor_, sstr);
  feasst_serialize_fstobj(particles_, sstr);
}

Cells::Cells(std::istream& sstr) {
  feasst_deserialize_version(sstr);
  feasst_deserialize(&label_, sstr);
  feasst_deserialize(&num_, sstr);
  feasst_deserialize(&neighbor_, sstr);
  feasst_deserialize_fstobj(&particles_, sstr);
}

}  // namespace feasst
