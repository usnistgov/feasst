
#ifndef FEASST_PATCH_PATCH_ANGLE_H_
#define FEASST_PATCH_PATCH_ANGLE_H_

#include "system/include/visit_model.h"
#include "math/include/utils_math.h"

namespace feasst {

class PatchAngle : public ModelParam {
 public:
  PatchAngle() { set_name("patch_angle"); }
};

class CosPatchAngle : public ModelParam {
 public:
  CosPatchAngle() { set_name("cos_patch_angle"); }
  CosPatchAngle(std::istream& istr) : ModelParam(istr) {}

  double compute(const int type, const ModelParams& model_params) override {
    const double angle = model_params.select("patch_angle")->value(type);
    return cos(degrees_to_radians(angle));
  }
};

}  // namespace feasst

#endif  // FEASST_PATCH_PATCH_ANGLE_H_
