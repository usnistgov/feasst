#include "utils/include/serialize.h"
#include "confinement/include/always_reject.h"

namespace feasst {

AlwaysReject::AlwaysReject(const argtype &args) : Criteria(args) {
  class_name_ = "AlwaysReject";
}

AlwaysReject::AlwaysReject(std::shared_ptr<Constraint> constraint,
    const argtype& args) : AlwaysReject(args) {
  add(constraint);
}

class MapAlwaysReject {
 public:
  MapAlwaysReject() {
    AlwaysReject().deserialize_map()["AlwaysReject"] = MakeAlwaysReject();
  }
};

static MapAlwaysReject mapper_ = MapAlwaysReject();

AlwaysReject::AlwaysReject(std::istream& istr) : Criteria(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 4204, "version mismatch: " << version);
}

void AlwaysReject::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_criteria_(ostr);
  feasst_serialize_version(4204, ostr);
}

}  // namespace feasst
