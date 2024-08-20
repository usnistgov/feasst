#include "utils/include/serialize.h"
#include "aniso/include/anisotropic.h"

namespace feasst {

class MapAnisotropic {
 public:
  MapAnisotropic() {
    auto obj = std::make_shared<Anisotropic>();
    obj->deserialize_map()["anisotropic"] = obj;
  }
};

static MapAnisotropic mapper_anisotropic_ = MapAnisotropic();

void Anisotropic::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(9568, ostr);
}

Anisotropic::Anisotropic(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9568, "mismatch version: " << version);
}

}  // namespace feasst
