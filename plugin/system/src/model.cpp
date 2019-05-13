#include "system/include/model.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Model> >& Model::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Model> >* ans =
     new std::map<std::string, std::shared_ptr<Model> >();
  return *ans;
}

}  // namespace feasst
