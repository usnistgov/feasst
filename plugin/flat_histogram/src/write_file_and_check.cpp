#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/file.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/write_file_and_check.h"

namespace feasst {

WriteFileAndCheck::WriteFileAndCheck(argtype * args) {
  class_name_ = "WriteFileAndCheck";
  sim_ = integer("sim", args);
  sim_start_ = integer("sim_start", args);
  sim_end_ = integer("sim_end", args);
  file_prefix_ = str("file_prefix", args);
  file_suffix_ = str("file_suffix", args);
  output_file_ = str("output_file", args);
}
WriteFileAndCheck::WriteFileAndCheck(argtype args) : WriteFileAndCheck(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(WriteFileAndCheck, argtype({{"sim", "0"}, {"sim_start", "0"},
  {"sim_end", "1"}, {"file_prefix", "a"}, {"file_suffix", "a"},
  {"output_file", "a"}}));

WriteFileAndCheck::WriteFileAndCheck(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 5744 && version <= 5744, "mismatch version: " << version);
  feasst_deserialize(&sim_, istr);
  feasst_deserialize(&sim_start_, istr);
  feasst_deserialize(&sim_end_, istr);
  feasst_deserialize(&file_prefix_, istr);
  feasst_deserialize(&file_suffix_, istr);
  feasst_deserialize(&output_file_, istr);
}

void WriteFileAndCheck::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(5744, ostr);
  feasst_serialize(sim_, ostr);
  feasst_serialize(sim_start_, ostr);
  feasst_serialize(sim_end_, ostr);
  feasst_serialize(file_prefix_, ostr);
  feasst_serialize(file_suffix_, ostr);
  feasst_serialize(output_file_, ostr);
}

void WriteFileAndCheck::run(MonteCarlo * mc) {
  std::ofstream file(file_prefix_ + str(sim_) + file_suffix_);
  bool all_exist = true;
  for (int index = sim_start_; index <= sim_end_; ++index) {
    if (!file_exists(file_prefix_ + str(index) + file_suffix_)) {
      all_exist = false;
      DEBUG("does not exist " << index);
    }
  }
  DEBUG("all_exist " << all_exist);
  if (all_exist) {
    std::ofstream file(output_file_);
  }
}

}  // namespace feasst
