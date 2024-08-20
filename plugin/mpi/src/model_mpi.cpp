#include <cmath>
#include "mpi.h"
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "mpi/include/thread_mpi.h"
#include "mpi/include/model_mpi.h"

namespace feasst {

class MapModelMPI {
 public:
  MapModelMPI() {
    ModelMPI().deserialize_map()["ModelMPI"] = MakeModelMPI();
  }
};

static MapModelMPI mapper_ = MapModelMPI();

ModelMPI::ModelMPI(argtype * args) {
  class_name_ = "ModelMPI";
  const double thres = dble("hard_sphere_threshold", args, 0.2);
  hard_sphere_threshold_sq_ = thres*thres;
}
ModelMPI::ModelMPI(argtype args) : ModelMPI(&args) {
  feasst_check_all_used(args);
}
ModelMPI::~ModelMPI() {
  if (thread_) {
    // Send termination signal as a squared distance less than 1
    double squared_distance = -1;
    MPI_Send(&squared_distance, 1 - thread_->thread(), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
  }
}

void ModelMPI::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_mpi_(ostr);
}

void ModelMPI::serialize_model_mpi_(std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(1056, ostr);
  feasst_serialize(hard_sphere_threshold_sq_, ostr);
}

ModelMPI::ModelMPI(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(1056 == version, version);
  feasst_deserialize(&hard_sphere_threshold_sq_, istr);
}

double ModelMPI::hard_sphere_threshold() const {
  return std::sqrt(hard_sphere_threshold_sq_);
}

//void ModelMPI::precompute(const ModelParams& existing) {
//  Model::precompute(existing);
//  int num = 1;
//  int initialized;
//  MPI_Initialized(&initialized);
//  if (!initialized) {
//    MPI_Init(NULL, NULL);
//  }
//  MPI_Comm_rank(MPI_COMM_WORLD, &thread_);
//  MPI_Comm_size(MPI_COMM_WORLD, &num);
//  ASSERT(thread_ == 0 || thread_ == 1, "thread:" << thread_ << " should be 0 or 1.");
//  ASSERT(num == 2, "num:" << num << " should be 2");
//}

double ModelMPI::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("squared_distance " << squared_distance);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  ASSERT(sigma_index() != -1, "err");
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  TRACE("sigma " << sigma);
  const double sigma_squared = sigma*sigma;
  if (squared_distance == 0 ||
      squared_distance < hard_sphere_threshold_sq_*sigma_squared) {
    TRACE("near inf");
    return NEAR_INFINITY;
  }
  if (!thread_) {
    thread_ = std::make_shared<ThreadMPI>();
    ASSERT(thread_->num() == 2, "num:" << thread_->num() << " should be 2");
  }
  MPI_Send(&squared_distance, 1 - thread_->thread(), MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
  MPI_Send(&type1, 1, MPI_INT, 1 - thread_->thread(), 0, MPI_COMM_WORLD);
  MPI_Send(&type2, 1, MPI_INT, 1 - thread_->thread(), 0, MPI_COMM_WORLD);
  MPI_Recv(&en_, 1, MPI_DOUBLE, 1 - thread_->thread(), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  return en_;
}

}  // namespace feasst
