#include "utils/include/serialize.h"
#include "steppers/include/read_config_from_file.h"

namespace feasst {

ReadConfigFromFile::ReadConfigFromFile(argtype * args) : ModifyUpdateOnly(args) {}
ReadConfigFromFile::ReadConfigFromFile(argtype args) : ReadConfigFromFile(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapReadConfigFromFile {
 public:
  MapReadConfigFromFile() {
    ReadConfigFromFile().deserialize_map()["ReadConfigFromFile"] = MakeReadConfigFromFile();
  }
};

static MapReadConfigFromFile mapper_ = MapReadConfigFromFile();

void ReadConfigFromFile::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  file_.open(file_name());
  ASSERT(file_.good(), "cannot open " << file_name());
}

void ReadConfigFromFile::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  Configuration * config = system->get_configuration();
  xyz_.load_frame(file_, config);
  if (file_.peek() == EOF) {
    criteria->set_num_iterations_to_complete(0);
  }
}

void ReadConfigFromFile::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(6782, ostr);
}

ReadConfigFromFile::ReadConfigFromFile(std::istream& istr) : ModifyUpdateOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6782, "version mismatch:" << version);
}

}  // namespace feasst
