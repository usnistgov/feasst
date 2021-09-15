
#ifndef FEASST_CONFIGURATION_SITE_H_
#define FEASST_CONFIGURATION_SITE_H_

#include <vector>
#include <string>
#include "math/include/position.h"
#include "configuration/include/typed_entity.h"
#include "configuration/include/properties.h"

namespace feasst {

/**
  Sites are used for interaction potentials.
  They contain information associated with their positions and identity.
  This Site base class contains only those site properties which are unique
  to this particular site and not the same for every site of the same type.
 */
class Site : public PropertiedEntity,
             public TypedEntity {
 public:
  Site();

  /// Set the Position.
  void set_position(const Position& position) { position_ = position; }

  /// Add to the position.
  void add_position(const Position& position) { position_.add(position); }

  /// Return the Position.
  const Position& position() const { return position_; }

  /// Return the Position in a given dimension.
  const double position(const int dimension) const {
    return position_.coord(dimension); }

  /// Construct with a position.
  Site(const Position position) : Site() { set_position(position); }

  /// Displace the Position of the Site.
  void displace(const Position displacement) {
    add_position(displacement); }

  /// Return true if the site is a director.
  /// These are used, for example, by patchy models.
  bool is_director() const { return is_director_; }

  /// Add a property with the given name and value.
  /// If the name is director, is_director() will return true.
  void add_property(const std::string name, const double value) override;

  /// Return true if the site is physically present.
  /// These are used, for example, in particle regrowth.
  bool is_physical() const { return is_physical_; }

  /// Set as physical/nonphysical (default: physical).
  void set_physical(const bool physical = true) { is_physical_ = physical; }

  /// Add a cell.
  void add_cell(const int cell) { cells_.push_back(cell); }

  /// Set a cell.
  void set_cell(const int index, const int cell) { cells_[index] = cell; }

  /// Return a cell.
  int cell(const int index) const;

  /// Return the number of cells.
  int num_cells() const { return static_cast<int>(cells_.size()); }

  void serialize(std::ostream& ostr) const;
  explicit Site(std::istream& istr);
  virtual ~Site() {}

 private:
  Position position_;
  bool is_director_ = false;
  bool is_physical_;
  std::vector<int> cells_;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_SITE_H_
