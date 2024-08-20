#include <cmath>
#include "mpi.h"
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "mpi/include/thread_mpi.h"
#include "mpi/include/mpi_place_holder.h"

namespace feasst {

class MapMPIPlaceHolder {
 public:
  MapMPIPlaceHolder() {
    MPIPlaceHolder().deserialize_map()["MPIPlaceHolder"] = MakeMPIPlaceHolder();
  }
};

static MapMPIPlaceHolder mapper_ = MapMPIPlaceHolder();

MPIPlaceHolder::MPIPlaceHolder(argtype * args) {
  class_name_ = "MPIPlaceHolder";
}
MPIPlaceHolder::MPIPlaceHolder(argtype args) : MPIPlaceHolder(&args) {
  feasst_check_all_used(args);
}

void MPIPlaceHolder::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(3257, ostr);
}

MPIPlaceHolder::MPIPlaceHolder(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(3257 == version, version);
}

double MPIPlaceHolder::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  FATAL("not implemented");
  return 0.;
}

}  // namespace feasst
