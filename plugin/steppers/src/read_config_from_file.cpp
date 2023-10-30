#include "utils/include/serialize.h"
#include "steppers/include/read_config_from_file.h"

namespace feasst {

ReadConfigFromFile::ReadConfigFromFile(argtype * args) : ModifyUpdateOnly(args) {
  input_file_ = str("input_file", args);
}
ReadConfigFromFile::ReadConfigFromFile(argtype args) : ReadConfigFromFile(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapReadConfigFromFile {
 public:
  MapReadConfigFromFile() {
    auto obj = MakeReadConfigFromFile({{"input_file", "placeholder"}});
    obj->deserialize_map()["ReadConfigFromFile"] = obj;
  }
};

static MapReadConfigFromFile mapper_ = MapReadConfigFromFile();

void ReadConfigFromFile::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  file_.open(input_file_);
  ASSERT(file_.good(), "cannot open " << input_file_);
}

void ReadConfigFromFile::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  INFO("ReadConfigFromFile::update");
  Configuration * config = system->get_configuration();
  if (xyz_.load_frame(file_, config)) {
    Acceptance acc_;
    criteria->update_state(*system, acc_);
    INFO("state " << criteria->state());
    if (file_.peek() == EOF) {
      INFO("setting complete");
      criteria->set_complete();
    }
  }
}

void ReadConfigFromFile::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(6782, ostr);
  feasst_serialize(input_file_, ostr);
}

ReadConfigFromFile::ReadConfigFromFile(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6782, "version mismatch:" << version);
  feasst_deserialize(&input_file_, istr);
}

}  // namespace feasst
