
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/write_model_params.h"

namespace feasst {

WriteModelParams::WriteModelParams(argtype * args) {
  class_name_ = "WriteModelParams";
  output_file_ = str("output_file", args);
  if (used("file_name", *args)) {
    WARN("WriteModelParams::file_name was renamed to output_file.");
    output_file_ = str("file_name", args);
  }
  potential_index_ = integer("potential_index", args, -1);
  reference_index_ = integer("reference_index", args, -1);
}
WriteModelParams::WriteModelParams(argtype args) : WriteModelParams(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapWriteModelParams {
 public:
  MapWriteModelParams() {
    auto obj = MakeWriteModelParams({{"output_file", "place_holder"}});
    obj->deserialize_map()["WriteModelParams"] = obj;
  }
};

static MapWriteModelParams mapper_WriteModelParams = MapWriteModelParams();

WriteModelParams::WriteModelParams(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2890, "mismatch version: " << version);
  feasst_deserialize(&output_file_, istr);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&reference_index_, istr);
}

void WriteModelParams::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2890, ostr);
  feasst_serialize(output_file_, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(reference_index_, ostr);
}

void WriteModelParams::run(MonteCarlo * mc) {
  std::ofstream file(output_file_);
  const Configuration& config = mc->configuration();
  if (potential_index_ == -1) {
    file << config.model_params().str();
  } else {
    if (reference_index_ == -1) {
      const Potential& potential = mc->system().potential(potential_index_);
      file << potential.model_params(config).str();
    } else {
      const Potential& potential = mc->system().reference(reference_index_,
                                                          potential_index_);
      file << potential.model_params(config).str();
    }
  }
}

}  // namespace feasst
