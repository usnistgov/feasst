#include "core/include/criteria_writer.h"

namespace feasst {

class MapCriteriaWriter {
 public:
  MapCriteriaWriter() {
    CriteriaWriter().deserialize_map()["CriteriaWriter"] = MakeCriteriaWriter();
  }
};

static MapCriteriaWriter mapper_ = MapCriteriaWriter();

}  // namespace feasst
