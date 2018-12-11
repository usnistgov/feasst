
#include <cmath>
#include "core/include/cells.h"
#include "core/include/utils_math.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"

namespace feasst {

//Cells::Cells() {
//  set_type();
//  set_group();
//}

void Cells::create(const double min_length,
                   const std::vector<double> side_lengths) {
  ASSERT(min_length > 1e-15, "min_length(" << min_length << ") too small");
  clear();
  for (double side_length : side_lengths) {
    num_.push_back(static_cast<int>(side_length/min_length));
//    lengths_.push_back(side_length/
//                            static_cast<double>(num_.back()));
  }
  if (num_total() <= std::pow(3, side_lengths.size())) {
    clear();
    return;
  }
  ASSERT(num_total() < 1e8, "too many cells");
  neighbor_.clear();
  neighbor_.resize(num_total());
  if (static_cast<int>(side_lengths.size()) == 3) {
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
  } else if (static_cast<int>(side_lengths.size()) == 2) {
    for (int xcell1 = 0; xcell1 < num(0); ++xcell1) {
    for (int ycell1 = 0; ycell1 < num(1); ++ycell1) {
      const int cell = id_({xcell1, ycell1});
      for (int xcell2 = xcell1 - 1; xcell2 <= xcell1 + 1; ++xcell2) {
      for (int ycell2 = ycell1 - 1; ycell2 <= ycell1 + 1; ++ycell2) {
        neighbor_[cell].push_back(id_({xcell2, ycell2}));
      }}
    }}
  } else {
    ASSERT(false, "unrecognized dimension(" << side_lengths.size() << ")");
  }
}

int Cells::num_total() const { return product(num_); }

//double Cells::length(const int dimension) const {
//  return lengths_[dimension];
//}
//
//double Cells::min_length() const { return minimum(lengths_); }

void Cells::clear() {
  num_.clear();
//  lengths_.clear();
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
  //for (const double& coord : scaled_coord) {
  for (int dim = 0; dim < static_cast<int>(cells.size()); ++dim) {
    ASSERT(std::abs(scaled_coord[dim]) < 0.5, "not scaled coordinates");
    cells[dim] = static_cast<int>(
      num_[dim]*(scaled_coord[dim] + 0.5)) % num_[dim];
    double prod = cells[dim];
    for (int dim2 = 0; dim2 < dim; ++dim2) {
      prod *= num_[dim2];
    }
    cell += prod;
  }
  ASSERT( (cell >= 0) && (cell < num_total()),
    "cell(" << cell << ")");
  return cell;
}

bool Cells::enabled() const {
  if (num_total() == 0) {
    return false;
  } else {
    return true;
  }
}

}  // namespace feasst
