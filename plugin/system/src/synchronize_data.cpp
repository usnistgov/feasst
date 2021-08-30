#include "utils/include/serialize.h"
#include "system/include/synchronize_data.h"

namespace feasst {

void SynchronizeData::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2462, ostr);
  feasst_serialize(dble_1D_, ostr);
  feasst_serialize(int_1D_, ostr);
  feasst_serialize(int64_1D_, ostr);
  feasst_serialize(dble_2D_, ostr);
  feasst_serialize(dble_3D_, ostr);
  feasst_serialize(dble_5D_, ostr);
  feasst_serialize(dble_6D_, ostr);
  feasst_serialize(vpvpvpvpv_, ostr);
  feasst_serialize(vvvpvpv_, ostr);
}

SynchronizeData::SynchronizeData(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2462, "mismatch version: " << version);
  feasst_deserialize(&dble_1D_, istr);
  feasst_deserialize(&int_1D_, istr);
  feasst_deserialize(&int64_1D_, istr);
  feasst_deserialize(&dble_2D_, istr);
  feasst_deserialize(&dble_3D_, istr);
  feasst_deserialize(&dble_5D_, istr);
  feasst_deserialize(&dble_6D_, istr);
  feasst_deserialize(&vpvpvpvpv_, istr);
  feasst_deserialize(&vvvpvpv_, istr);
}

}  // namespace feasst
