#include "core/include/site.h"
#include "core/include/debug.h"

namespace feasst {

void Site::add_property(const std::string name, const double value) {
  PropertiedEntity::add_property(name, value);
  if (name == "director") {
    is_director_ = true;
  }
}

}  // namespace feasst
