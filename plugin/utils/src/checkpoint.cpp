#include "utils/include/checkpoint.h"
#include "utils/include/utils_io.h"
#include "utils/include/serialize.h"

namespace feasst {

Checkpoint::Checkpoint(const argtype &args) {
  args_.init(args);
  num_hours_ = args_.key("num_hours").dflt("12").dble();
  file_name_ = args_.key("file_name").dflt(" ").str();
}

void Checkpoint::serialize(std::ostream& ostr) const {
  feasst_serialize_version(223, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(previous_hours_, ostr);
  feasst_serialize(num_hours_, ostr);
}

Checkpoint::Checkpoint(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 223, "version mismatch: " << version);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&previous_hours_, istr);
  feasst_deserialize(&num_hours_, istr);
}

}  // namespace feasst
