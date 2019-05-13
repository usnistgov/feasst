#include "system/include/long_range_corrections.h"

namespace feasst {

class MapLongRangeCorrections {
 public:
  MapLongRangeCorrections() {
    LongRangeCorrections().deserialize_map()["LongRangeCorrections"] =
      MakeLongRangeCorrections();
  }
};

static MapLongRangeCorrections mapper_ = MapLongRangeCorrections();

}  // namespace feasst
