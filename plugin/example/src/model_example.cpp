#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/model_params.h"
#include "example/include/model_example.h"

namespace feasst {

class MapModelExample {
 public:
  MapModelExample() {
    auto obj = MakeModelExample({{"example_argument", "1"}});
    obj->deserialize_map()["ModelExample"] = obj;
  }
};

static MapModelExample mapper_ = MapModelExample();

/*
  Here input file arguments are parsed and stored as private member variables.

  The "example_argument" is of type double precision, and its default value is 1.
  For a different type, see avaible options in feasst/plugin/utils/include/arguments.h
  These include str, dble, flt, integer and boolean.

  If no default argument is provided, then the argument must always be provided,
  even in the MapModelExample() above.
 */
ModelExample::ModelExample(argtype * args) {
  class_name_ = "ModelExample";
  example_argument_ = dble("example_argument", args, 1);
}
ModelExample::ModelExample(argtype args) : ModelExample(&args) {
  FEASST_CHECK_ALL_USED(args);
}

/*
  Here the energy is finally calculated between two sites.
  As an example, lets implement the Jagla potential which has a hard sphere
  at sigma, and then a linear ramp from sigma to cutoff starting with an
  energy of epsilon at sigma and ending with an energy of zero at cutoff.
 */
double ModelExample::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  INFO("here");
  // TRACE only prints at a VERBOSE_LEVEL = 0 in /feasst/plugin/utils/include/debug.h
  // Replace TRACE with INFO to print at the default VERBOSE_LEVEL of 3.
  TRACE("squared_distance " << squared_distance);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
  TRACE("sigma " << sigma);
  const double sigma_squared = sigma*sigma;
  if (squared_distance < sigma_squared) {
    return NEAR_INFINITY;
  }
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  const double distance = std::sqrt(squared_distance);
  const double z = (cutoff - distance)/(cutoff - sigma);
  const double epsilon = model_params.select(epsilon_index()).mixed_values()[type1][type2];
  const double en = epsilon*z;
  TRACE("en " << en);
  return en;
}

/**
  If you want your custom model to be serializable, you must
  serialize any data member variables.
 */
void ModelExample::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(5023, ostr);
  feasst_serialize(example_argument_, ostr);
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
  feasst_deserialize(&example_argument_, istr);
}

}  // namespace feasst
