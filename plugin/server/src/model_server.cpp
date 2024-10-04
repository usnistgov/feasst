#include <cmath>
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "server/include/server.h"
#include "server/include/model_server.h"

namespace feasst {

FEASST_MAPPER(ModelServer,);

ModelServer::ModelServer(argtype * args) {
  class_name_ = "ModelServer";
  const double thres = dble("hard_sphere_threshold", args, 0.2);
  hard_sphere_threshold_sq_ = thres*thres;
  server_ = std::make_unique<Server>(args);
}
ModelServer::ModelServer(argtype args) : ModelServer(&args) {
  feasst_check_all_used(args);
}
ModelServer::~ModelServer() {}

void ModelServer::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_server_(ostr);
}

void ModelServer::serialize_model_server_(std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(2367, ostr);
  feasst_serialize(hard_sphere_threshold_sq_, ostr);
  feasst_serialize(server_, ostr);
}

ModelServer::ModelServer(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2367 == version, version);
  feasst_deserialize(&hard_sphere_threshold_sq_, istr);
  feasst_deserialize(server_, istr);
}

double ModelServer::hard_sphere_threshold() const {
  return std::sqrt(hard_sphere_threshold_sq_);
}

double ModelServer::energy(
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
  ASSERT(server_, "error");
  if (!server_->bound()) {
    server_->bind_listen_accept();
  }
  std::stringstream ss;
  ss << squared_distance << "," << type1 << "," << type2;
  server_->send(ss.str());
  //sprintf(server_->get_buffer(), "%f", squared_distance);
  //server_->buffer_() = sprintf;
  //std::ftoa(squared_distance));
  //server_->send_buffer();
  const int size = server_->receive();
  TRACE(server_->buffer());
  ASSERT(size > 0, "error");
  return std::stod(server_->buffer());
}

}  // namespace feasst
