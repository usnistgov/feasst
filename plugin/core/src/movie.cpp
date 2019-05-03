#include "core/include/movie.h"

namespace feasst {

class MapMovie {
 public:
  MapMovie() {
    auto obj = MakeMovie({{"file_name", "place_holder"}});
    obj->deserialize_map()["Movie"] = obj;
  }
};

static MapMovie mapper_ = MapMovie();

}  // namespace feasst
