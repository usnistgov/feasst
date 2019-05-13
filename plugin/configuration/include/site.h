
#ifndef FEASST_CONFIGURATION_SITE_H_
#define FEASST_CONFIGURATION_SITE_H_

#include <vector>
#include <utility>
#include <string>
#include "configuration/include/typed_entity.h"
#include "math/include/position.h"
#include "utils/include/utils_io.h"
#include "configuration/include/properties.h"

namespace feasst {

/**
  Sites are used for interaction potentials.
  They contain information associated with their positions and identity.
  This Site base class contains only those site properties which are unique
  to this particular site and not the same for every site of the same type.
 */
class Site : public PropertiedEntity,
             public TypedEntity,
             public SpatialEntity {
 public:
  Site() {}

  /// Displace the Position of the Site.
  void displace(const Position displacement) {
    add_position(displacement);
  }

  /// Return true if the site is a director.
  /// These are used, for example, by patchy models.
  bool is_director() const { return is_director_; }

  /// Assign as director if a site has a property named director.
  void add_property(const std::string name, const double value) override;

  void serialize(std::ostream& ostr) const;
  Site(std::istream& istr);
  virtual ~Site() {}

 private:
  bool is_director_ = false;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_SITE_H_
