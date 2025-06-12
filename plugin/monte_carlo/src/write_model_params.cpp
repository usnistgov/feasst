#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "utils/include/timer.h"
#include "math/include/utils_math.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/potential.h"
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
  if (used("reference_index", *args)) {
    WARN("Deprecated WriteModelParams::reference_index->ref.");
  }
  reference_index_ = integer("reference_index", args, -1);
  ref_ = str("ref", args, "");
  config_ = str("config", args, "0");
}
WriteModelParams::WriteModelParams(argtype args) : WriteModelParams(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(WriteModelParams, argtype({{"output_file", "place_holder"}}));

WriteModelParams::WriteModelParams(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 2890 && version <= 2891, "mismatch version: " << version);
  feasst_deserialize(&output_file_, istr);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&reference_index_, istr);
  if (version >= 2891) {
    feasst_deserialize(&ref_, istr);
    feasst_deserialize(&config_, istr);
  }
}

void WriteModelParams::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2891, ostr);
  feasst_serialize(output_file_, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(reference_index_, ostr);
  feasst_serialize(ref_, ostr);
  feasst_serialize(config_, ostr);
}

void WriteModelParams::run(MonteCarlo * mc) {
  std::ofstream file(output_file_);
  const int config_index = mc->system().configuration_index(config_);
  const Configuration& config = mc->system().configuration(config_);
  std::vector<std::string> site_type_names;
  config.site_type_names(&site_type_names);
  if (potential_index_ == -1) {
    file << config.model_params().str(&site_type_names);
  } else {
    if (!ref_.empty()) {
      reference_index_ = mc->system().reference_index(config_index, ref_);
    }
    if (reference_index_ == -1) {
      const Potential& potential = mc->system().potential(potential_index_);
      file << potential.model_params(config).str(&site_type_names);
    } else {
      const Potential& potential = mc->system().reference(reference_index_,
                                                          potential_index_);
      file << potential.model_params(config).str(&site_type_names);
    }
  }
}

}  // namespace feasst
