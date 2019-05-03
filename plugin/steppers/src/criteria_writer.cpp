#include "steppers/include/criteria_writer.h"

namespace feasst {

CriteriaWriter::CriteriaWriter(const argtype &args) : AnalyzeWriteOnly(args) {}

class MapCriteriaWriter {
 public:
  MapCriteriaWriter() {
    CriteriaWriter().deserialize_map()["CriteriaWriter"] = MakeCriteriaWriter();
  }
};

static MapCriteriaWriter mapper_ = MapCriteriaWriter();

}  // namespace feasst
