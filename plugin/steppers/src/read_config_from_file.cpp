#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/read_config_from_file.h"

namespace feasst {

ReadConfigFromFile::ReadConfigFromFile(argtype * args) : ModifyUpdateOnly(args) {
  input_file_ = str("input_file", args);
  xyz_ = FileXYZ(args);
}
ReadConfigFromFile::ReadConfigFromFile(argtype args) : ReadConfigFromFile(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(ReadConfigFromFile, argtype({{"input_file", "placeholder"}}));

void ReadConfigFromFile::load_(Criteria * criteria, System * system) {
  if (set_complete_next_update_) {
    DEBUG("setting complete");
    criteria->set_complete();
    return;
  }
  Configuration * config = system->get_configuration();
  if (xyz_.load_frame(file_, config)) {
    Acceptance acc_;
    criteria->update_state(*system, acc_);
    DEBUG("state " << criteria->state());
    if (file_.peek() == EOF) {
      set_complete_next_update_ = true;
    }
  }
}

void ReadConfigFromFile::initialize(MonteCarlo * mc) {
  Modify::initialize(mc);
  file_.open(input_file_);
  ASSERT(file_.good(), "cannot open " << input_file_);
  load_(mc->get_criteria(), mc->get_system());
}

void ReadConfigFromFile::update(MonteCarlo * mc) {
  DEBUG("ReadConfigFromFile::update");
  load_(mc->get_criteria(), mc->get_system());
}

void ReadConfigFromFile::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(6782, ostr);
  feasst_serialize(input_file_, ostr);
  feasst_serialize(set_complete_next_update_, ostr);
  feasst_serialize_fstobj(xyz_, ostr);
}

ReadConfigFromFile::ReadConfigFromFile(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6782, "version mismatch:" << version);
  feasst_deserialize(&input_file_, istr);
  feasst_deserialize(&set_complete_next_update_, istr);
  feasst_deserialize_fstobj(&xyz_, istr);
}

}  // namespace feasst
