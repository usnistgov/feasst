#include "core/include/movie.h"

namespace feasst {

class MapMovie {
 public:
  MapMovie() {
    Movie().deserialize_map()["Movie"] = MakeMovie();
  }
};

static MapMovie mapper_ = MapMovie();

}  // namespace feasst
