
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
  double compute(const int type, const ModelParams& model_params) override;
  CosPatchAngle(std::istream& istr) : ModelParam(istr) {}
};

class Director : public ModelParam {
 public:
  Director() { set_name("director"); }
  double compute(const int type, const ModelParams& model_params) override;
  Director(std::istream& istr) : ModelParam(istr) {}
};

}  // namespace feasst

#endif  // FEASST_PATCH_PATCH_ANGLE_H_
