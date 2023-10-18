
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/write_checkpoint.h"

namespace feasst {

WriteCheckpoint::WriteCheckpoint(argtype * args) {
  class_name_ = "WriteCheckpoint";
}
WriteCheckpoint::WriteCheckpoint(argtype args) : WriteCheckpoint(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapWriteCheckpoint {
 public:
  MapWriteCheckpoint() {
    auto obj = MakeWriteCheckpoint();
    obj->deserialize_map()["WriteCheckpoint"] = obj;
  }
};

static MapWriteCheckpoint mapper_WriteCheckpoint = MapWriteCheckpoint();

WriteCheckpoint::WriteCheckpoint(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5694, "mismatch version: " << version);
}

void WriteCheckpoint::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(5694, ostr);
}

void WriteCheckpoint::run(MonteCarlo * mc) {
  mc->write_checkpoint();
}

}  // namespace feasst
