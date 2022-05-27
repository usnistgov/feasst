#include "steppers/include/criteria_writer.h"
#include "utils/include/serialize.h"

namespace feasst {

CriteriaWriter::CriteriaWriter(argtype * args) : AnalyzeWriteOnly(args) {}
CriteriaWriter::CriteriaWriter(argtype args) : CriteriaWriter(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapCriteriaWriter {
 public:
  MapCriteriaWriter() {
    CriteriaWriter().deserialize_map()["CriteriaWriter"] = MakeCriteriaWriter();
  }
};

static MapCriteriaWriter mapper_ = MapCriteriaWriter();

std::string CriteriaWriter::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  // ensure the following order matches the header from initialization.
  std::stringstream ss;
  ss << criteria.write() << std::endl;
  return ss.str();
}

void CriteriaWriter::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(993, ostr);
}

CriteriaWriter::CriteriaWriter(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 993, "version mismatch:" << version);
}

}  // namespace feasst
