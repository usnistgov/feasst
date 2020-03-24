
#ifndef FEASST_PATCH_PATCH_ANGLE_H_
#define FEASST_PATCH_PATCH_ANGLE_H_

#include "configuration/include/model_params.h"

namespace feasst {

class PatchAngle : public ModelParam {
 public:
  PatchAngle() { set_name("patch_angle"); }
};

class CosPatchAngle : public ModelParam {
 public:
  CosPatchAngle() { set_name("cos_patch_angle"); }
  CosPatchAngle(std::istream& istr) : ModelParam(istr) {}

  double compute(const int type, const ModelParams& model_params) override;
};

}  // namespace feasst

#endif  // FEASST_PATCH_PATCH_ANGLE_H_
