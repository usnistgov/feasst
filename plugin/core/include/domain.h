
#ifndef FEASST_CORE_DOMAIN_H_
#define FEASST_CORE_DOMAIN_H_

#include <vector>
#include "core/include/position.h"
#include "core/include/random.h"
#include "core/include/cells.h"

namespace feasst {

/**
  A Domain represents the spatial boundaries and constraints imposed upon the
  positions of the particles.

  By convention, the origin is located at the center of the domain.

  By default, periodicity in each dimension is enabled.
  But it may be disabled manually.

  HWH Note: Cells are specific to cuboid domain.
*/
class Domain {
 public:
  /// Get the side lengths.
  Position side_length() const { return side_length_; }

  /// Get the side length.
  double side_length(const int dimension) const {
    return side_length_.coord(dimension);
  }

  /// Set the side lengths.
  void set_side_length(const Position& side_length) {
    side_length_ = side_length;
    periodic_.resize(side_length_.size(), true);
  }

  /// Set the side length.
  void set_side_length(const int dimension, const double length) {
    side_length_.set_coord(dimension, length);
  }

  /// Disable periodicity in a particular dimension.
  void disable(const int dimension) { periodic_[dimension] = false; }

  /// Return the dimensionality.
  int dimension() const { return side_length_.size(); }

  /// Return the volume.
  double volume() const;

  /// Return the shift necessary to wrap the position.
  virtual Position shift(const Position& position) const = 0;

  /// Wrap the input position.
  void wrap(Position * position) const;

  /// Return random position within boundaries.
  /// Require random number generator to keep method constant and prevent
  /// seeding to the same value over and over again.
  virtual Position random_position(Random * random) const = 0;

  /// Return the minimal side length
  double min_side_length() const;

  /// HWH implement check_size

  /// Initialize the cells according to the minimum side length.
  void init_cells(const double min_length);

  /// Return the cells.
  const std::vector<Cells>& cells() const { return cells_; }

  /// Return the unique number cell for the position.
  int cell_id(const Position& position, const Cells& cells) const;

  virtual ~Domain() {}

 protected:
  Position side_length_;
  std::vector<bool> periodic_;

  /// Cell lists
  std::vector<Cells> cells_;
};

}  // namespace feasst

#endif  // FEASST_CORE_DOMAIN_H_
