#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/model_params.h"
#include "example/include/model_example.h"

namespace feasst {

FEASST_MAPPER(ModelExample, argtype({{"num_discretized_steps", "1"}}));

/*
  In this constructor, input file arguments are parsed and stored as private
  member variables.

  The "num_discretized_steps" is of type integer, and its default value is 0.
  For a different type, see avaible options in feasst/plugin/utils/include/arguments.h
  These include str, dble, flt, integer and boolean.

  If no default argument is provided, then the argument must always be provided,
  even in the MapModelExample() above.
 */
ModelExample::ModelExample(argtype * args) {
  class_name_ = "ModelExample";
  num_discretized_steps_ = integer(
    "num_discretized_steps",  // name of argument
    args,                     // args pointer
    0);                       // (optional) default value if no argument given.
}
ModelExample::ModelExample(argtype args) : ModelExample(&args) {
  feasst_check_all_used(args);
}

void ModelExample::precompute(const ModelParams& existing) {
  Model::precompute(existing);  // find sigma, epsilon, cutoff, charge index
  lambda_index_ = existing.index("lambda");
  gamma_index_ = existing.index("gamma");
}

/*
  The energy is calculated between two sites.
 */
double ModelExample::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  // TRACE prints at a VERBOSE_LEVEL = 0 in /feasst/plugin/utils/include/debug.h
  // Replace TRACE with INFO to print at the default VERBOSE_LEVEL of 3.
  TRACE("squared_distance " << squared_distance);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  ASSERT(num_discretized_steps_ == 0,
    "As an exercise, consider implementing num_discretized_steps to compare " <<
    "with the results of discrete molecular dynamics simulations.");
  const double sigma =
    model_params.select(sigma_index()).mixed_value(type1, type2);
  TRACE("sigma " << sigma);
  if (squared_distance < sigma*sigma) {
    return NEAR_INFINITY;
  }
  const double epsilon =
    model_params.select(epsilon_index()).mixed_value(type1, type2);
  const double cutoff =
    model_params.select(cutoff_index()).mixed_value(type1, type2);
  const double distance = std::sqrt(squared_distance);
  const double lambda =
    model_params.select(lambda_index_).mixed_value(type1, type2);
  if (squared_distance < lambda*lambda) {
    const double gamma =
      model_params.select(gamma_index_).mixed_value(type1, type2);
    const double en = (gamma*(lambda - distance) - epsilon*(distance - sigma))/
                      (lambda - sigma);
    TRACE("en " << en);
    return en;
  } else {
    const double en = -epsilon*(cutoff - distance)/(cutoff - lambda);
    TRACE("en " << en);
    return en;
  }
  // Note that VisitModel should never call energy beyond its cutoff distance,
  // so for optimization, there is no need for the check here.
}

/**
  If you want your custom model to be serializable, you must
  serialize any data member variables.
 */
void ModelExample::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(5023, ostr);
  feasst_serialize(num_discretized_steps_, ostr);
  feasst_serialize(lambda_index_, ostr);
  feasst_serialize(gamma_index_, ostr);
}

/*
  The same data member variables must be deserialized in the exact same
  order as they were serialized above.
  Also, consider changing the version number '5023' in both cases to some
  other random four digit number.
  This helps test the serialization string for errors and add backwards
  compatibility for reading serializations of older versions.
 */
ModelExample::ModelExample(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5023, "unrecognized version: " << version);
  feasst_deserialize(&num_discretized_steps_, istr);
  feasst_deserialize(&lambda_index_, istr);
  feasst_deserialize(&gamma_index_, istr);
}

}  // namespace feasst
