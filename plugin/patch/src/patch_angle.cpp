#include <cmath>
#include "patch/include/patch_angle.h"
#include "math/include/utils_math.h"

namespace feasst {

double CosPatchAngle::compute(const int type, const ModelParams& model_params) {
  const double angle = model_params.select("patch_angle")->value(type);
  return std::cos(degrees_to_radians(angle));
}

}  // namespace feasst
