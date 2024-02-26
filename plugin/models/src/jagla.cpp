#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/model_params.h"
#include "example/include/model_param_example.h"
#include "models/include/jagla.h"

namespace feasst {

class MapJagla {
 public:
  MapJagla() {
    auto obj = MakeJagla({{"num_discretized_steps", "0"}});
    obj->deserialize_map()["Jagla"] = obj;
  }
};

static MapJagla mapper_ = MapJagla();

Jagla::Jagla(argtype * args) {
  class_name_ = "Jagla";
  num_discretized_steps_ = integer("num_discretized_steps", args, 0);
  ASSERT(num_discretized_steps_ == 0, "num_discretized_steps: " <<
    num_discretized_steps_ << " not implemented");
}
Jagla::Jagla(argtype args) : Jagla(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void Jagla::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  lambda_index_ = existing.index("lambda");
  gamma_index_ = existing.index("gamma");
}

double Jagla::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("squared_distance " << squared_distance);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  TRACE("sigma " << sigma);
  if (squared_distance < sigma*sigma) {
    return NEAR_INFINITY;
  }
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  const double distance = std::sqrt(squared_distance);
  const double lambda = model_params.select(lambda_index_).mixed_values()[type1][type2];
  if (squared_distance < lambda*lambda) {
    const double gamma = model_params.select(gamma_index_).mixed_values()[type1][type2];
    const double en = (gamma*(lambda - distance) - epsilon*(distance - sigma))/(lambda - sigma);
    TRACE("en " << en);
    return en;
  } else {
    const double en = -epsilon*(cutoff - distance)/(cutoff - lambda);
    TRACE("en " << en);
    return en;
  }
}

void Jagla::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(6349, ostr);
  feasst_serialize(num_discretized_steps_, ostr);
  feasst_serialize(lambda_index_, ostr);
  feasst_serialize(gamma_index_, ostr);
}

Jagla::Jagla(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6349, "unrecognized version: " << version);
  feasst_deserialize(&num_discretized_steps_, istr);
  feasst_deserialize(&lambda_index_, istr);
  feasst_deserialize(&gamma_index_, istr);
}

}  // namespace feasst
