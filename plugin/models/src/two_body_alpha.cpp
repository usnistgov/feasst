#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "models/include/two_body_alpha.h"

namespace feasst {

TwoBodyAlpha::TwoBodyAlpha(argtype * args) {
  class_name_ = "TwoBodyAlpha";

  DEBUG("parse alpha and s");
  std::string start;
  std::stringstream key;
  start.assign("alpha");
  int index = 0;
  key << start << index;
  DEBUG("key " << key.str());
  while (used(key.str(), *args)) {
    alpha_.push_back(dble(key.str(), args));
    key.str("");
    key << "s" << index;
    DEBUG("key " << key.str());
    s_.push_back(dble(key.str(), args));
    ++index;
    key.str("");
    key << start << index;
    DEBUG("key " << key.str());
    DEBUG("args " << str(*args));
  }
}
TwoBodyAlpha::TwoBodyAlpha(argtype args) : TwoBodyAlpha(&args) {
  feasst_check_all_used(args);
}

FEASST_MAPPER(TwoBodyAlpha, argtype({{"alpha0", "1"}, {"s0", "1"}}));

void TwoBodyAlpha::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_two_body_alpha_(ostr);
}

void TwoBodyAlpha::serialize_two_body_alpha_(
    std::ostream& ostr) const {
  serialize_model_(ostr);
  feasst_serialize_version(3479, ostr);
  feasst_serialize(alpha_, ostr);
  feasst_serialize(s_, ostr);
}

TwoBodyAlpha::TwoBodyAlpha(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(3479 == version, version);
  feasst_deserialize(&alpha_, istr);
  feasst_deserialize(&s_, istr);
}

double TwoBodyAlpha::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  const double sigma_squared = sigma*sigma;
  TRACE("squared_distance " << squared_distance);
  TRACE("sigma_squared " << sigma_squared);
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double distance = std::sqrt(squared_distance);
  TRACE("distance " << distance);
  const double s_inv_dist = sigma/distance;
  double en = 0.;
  for (int index = 0; index < static_cast<int>(alpha_.size()); ++index) {
    const double alpha = alpha_[index];
    const double s = s_[index];
    en += s*epsilon*std::pow(s_inv_dist, alpha);
  }
  return en;
}

}  // namespace feasst
