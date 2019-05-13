
#include "utils/include/debug.h"
#include "monte_carlo/include/modify.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Modify> >& Modify::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Modify> >* ans =
     new std::map<std::string, std::shared_ptr<Modify> >();
  return *ans;
}

void Modify::serialize(std::ostream& ostr) const { ERROR("not implemented"); }

std::shared_ptr<Modify> Modify::create(std::istream& istr) const {
  ERROR("not implemented");
}

std::shared_ptr<Modify> Modify::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(),
    istr,
    // rewind reading of class name
    true);
}

}  // namespace feasst
