#include "system/include/ideal_gas.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapIdealGas {
 public:
  MapIdealGas() {
    IdealGas().deserialize_map()["IdealGas"] = MakeIdealGas();
  }
};

static MapIdealGas mapper_ = MapIdealGas();

void IdealGas::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_ideal_gas_(ostr);
}

void IdealGas::serialize_ideal_gas_(std::ostream& ostr) const {
  feasst_serialize_version(293, ostr);
}

IdealGas::IdealGas(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(293 == version, version);
}

}  // namespace feasst
