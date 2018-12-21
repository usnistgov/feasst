
#ifndef FEASST_CORE_SITE_H_
#define FEASST_CORE_SITE_H_

#include <vector>
#include <utility>
#include <string>
#include "core/include/typed_entity.h"
#include "core/include/position.h"
#include "core/include/utils_io.h"
#include "core/include/properties.h"

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
  /// Displace the Position of the Site.
  void displace(const Position displacement) {
    add_position(displacement);
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_SITE_H_
