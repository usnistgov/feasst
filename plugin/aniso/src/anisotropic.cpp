#include "utils/include/serialize.h"
#include "aniso/include/anisotropic.h"

namespace feasst {

FEASST_MAPPER_RENAME(Anisotropic, anisotropic,);

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
