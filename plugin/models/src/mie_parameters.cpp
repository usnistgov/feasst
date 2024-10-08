#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "models/include/mie_parameters.h"

namespace feasst {

double MieLambdaR::mix_(const double value1, const double value2) {
  return std::sqrt((value1 - 3.)*(value2 - 3.)) + 3.;
}

FEASST_MAPPER_RENAME(MieLambdaR, mie_lambda_r,);

void MieLambdaR::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(7083, ostr);
}

MieLambdaR::MieLambdaR(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7083, "mismatch version: " << version);
}

double MieLambdaA::mix_(const double value1, const double value2) {
  const double mix = std::sqrt((value1 - 3.)*(value2 - 3.)) + 3.;
  return mix;
}

FEASST_MAPPER_RENAME(MieLambdaA, mie_lambda_a,);

void MieLambdaA::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2670, ostr);
}

MieLambdaA::MieLambdaA(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2670, "mismatch version: " << version);
}

FEASST_MAPPER_RENAME(MiePrefactor, mie_prefactor,);

void MiePrefactor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(7964, ostr);
}

MiePrefactor::MiePrefactor(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7964, "mismatch version: " << version);
}

double MiePrefactor::compute(const int type1, const int type2, const ModelParams& model_params) {
  const double eps = model_params.select("epsilon").mixed_value(type1, type2);
  const double n = model_params.select("mie_lambda_r").mixed_value(type1, type2);
  const double m = model_params.select("mie_lambda_a").mixed_value(type1, type2);
  return eps*n/(n-m)*std::pow(n/m, m/(n-m));
}

//double MieIdealDeviation::mix_(const double value1, const double value2) {
//  return 0.;
//}
//
//FEASST_MAPPER_RENAME(MieIdealDeviation, mie_ideal_deviation,);
//
//void MieIdealDeviation::serialize(std::ostream& ostr) const {
//  ostr << class_name_ << " ";
//  serialize_model_param_(ostr);
//  feasst_serialize_version(2073, ostr);
//}
//
//MieIdealDeviation::MieIdealDeviation(std::istream& istr) : ModelParam(istr) {
//  const int version = feasst_deserialize_version(istr);
//  ASSERT(version == 2073, "mismatch version: " << version);
//}

}  // namespace feasst
