#include "utils/include/checkpoint.h"
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"

namespace feasst {

Checkpoint::Checkpoint(argtype args) {
  num_hours_ = dble("num_hours", &args, 1.);
  num_hours_terminate_ = dble("num_hours_terminate", &args, -1);
  checkpoint_file_ = str("checkpoint_file", &args, " ");
  if (used("file_name", args)) {
    WARN("Checkpoint::file_name renamed to checkpoint_file.");
    checkpoint_file_ = str("file_name", &args);
  }
  writes_per_backup_ = integer("writes_per_backup", &args, -1);
  first_hours_ = cpu_hours();
  feasst_check_all_used(args);
}

void Checkpoint::serialize(std::ostream& ostr) const {
  feasst_serialize_version(223, ostr);
  feasst_serialize(checkpoint_file_, ostr);
  feasst_serialize(num_hours_, ostr);
  feasst_serialize(num_hours_terminate_, ostr);
  feasst_serialize(writes_per_backup_, ostr);
  feasst_serialize(previous_backup_, ostr);
}

Checkpoint::Checkpoint(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 223, "version mismatch: " << version);
  feasst_deserialize(&checkpoint_file_, istr);
  feasst_deserialize(&num_hours_, istr);
  feasst_deserialize(&num_hours_terminate_, istr);
  feasst_deserialize(&writes_per_backup_, istr);
  feasst_deserialize(&previous_backup_, istr);
  first_hours_ = cpu_hours();
}

}  // namespace feasst
