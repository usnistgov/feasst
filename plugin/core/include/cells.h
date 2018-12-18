
#ifndef FEASST_CORE_CELLS_H_
#define FEASST_CORE_CELLS_H_

#include <string>
#include <vector>
#include "core/include/properties.h"

namespace feasst {

/**
  Divide a cuboid domain into cells.
 */
class Cells : public PropertiedEntity {
 public:
  /// Create the number, length and neighbors.
  /// By default, abort if there aren't more than \f$3^D\f$ cells,
  /// where D is the dimension.
  // HWH: better optimize method of building list of neighboring cells.
  // HWH: currently too slow to NPT, frequent rebuilds or large systems.
  // HWH: consider elongated boxes for minimal requirement.
  void create(const double min_length, const std::vector<double> side_lengths);

  /// Return the number.
  int num_total() const;

  /// Return if cells are enabled.
  bool enabled() const;

  /// Return the number in a dimension.
  int num(const int dimension) const { return num_[dimension]; }

  /// Return the number.
  std::vector<int> num() const { return num_; }

  /// Clear all private member data.
  void clear();

  /// Return the neighbors.
  /// The first index is the cell.
  /// The second is a list of neighboring cells (including self).
  const std::vector<std::vector<int> >& neighbor() const { return neighbor_; }

  /// Return the unique number cell in which the scaled coordinate resides.
  /// Scaled coordinates are positions divided by the respective domain size.
  int id(const std::vector<double>& scaled_coord) const;

  /// Return the label.
  std::string label() const { return label_; }

  /// Set the label.
  void set_label(const std::string label) { label_ = label; }

 private:
  std::vector<int> num_;
  std::vector<std::vector<int> > neighbor_;
  std::string label_;

  /// Return the unique cell id number for a given cell vector.
  int id_(std::vector<int> position);
};

}  // namespace feasst

#endif  // FEASST_CORE_CELLS_H_
