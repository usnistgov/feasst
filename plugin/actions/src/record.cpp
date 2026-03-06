#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "configuration/include/configuration.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/file_vmd.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "actions/include/record.h"

namespace feasst {

Record::Record(argtype * args) {
  class_name_ = "Record";
  config_ = str("config", args, "0");
  save_positions_ = str("save_positions", args, "");
  load_positions_ = str("load_positions", args, "");
}
Record::Record(argtype args) : Record(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(Record,);

Record::Record(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 8364 && version <= 8364, "mismatch version: " << version);
  feasst_deserialize(&config_, istr);
  feasst_deserialize(&save_positions_, istr);
  feasst_deserialize(&load_positions_, istr);
}

void Record::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(8364, ostr);
  feasst_serialize(config_, ostr);
  feasst_serialize(save_positions_, ostr);
  feasst_serialize(load_positions_, ostr);
}

void Record::run(MonteCarlo * mc) {
  if (!save_positions_.empty()) {
    const Configuration& config = mc->system().configuration(config_);
    FileXYZ xyz;
    xyz.write(save_positions_, config);
    FileVMD vmd;
    vmd.write(save_positions_ + ".vmd", config, save_positions_);
  }
  if (!load_positions_.empty()) {
    const int config_index = mc->system().configuration_index(config_);
    FileXYZ xyz;
    xyz.load(load_positions_, mc->get_system()->get_configuration(config_index));
    mc->initialize_criteria();
  }
}

}  // namespace feasst
