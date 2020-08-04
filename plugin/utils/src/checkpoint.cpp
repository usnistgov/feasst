#include "utils/include/checkpoint.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"

namespace feasst {

Checkpoint::Checkpoint(const argtype &args) {
  Arguments args_(args);
  num_hours_ = args_.key("num_hours").dflt("1.").dble();
  num_hours_terminate_ = args_.key("num_hours_terminate").dflt("-1").dble();
  file_name_ = args_.key("file_name").dflt(" ").str();
  first_hours_ = cpu_hours();
}

void Checkpoint::serialize(std::ostream& ostr) const {
  feasst_serialize_version(223, ostr);
  feasst_serialize(file_name_, ostr);
  feasst_serialize(num_hours_, ostr);
  feasst_serialize(num_hours_terminate_, ostr);
}

Checkpoint::Checkpoint(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 223, "version mismatch: " << version);
  feasst_deserialize(&file_name_, istr);
  feasst_deserialize(&num_hours_, istr);
  feasst_deserialize(&num_hours_terminate_, istr);
  first_hours_ = cpu_hours();
}

}  // namespace feasst
