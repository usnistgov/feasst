
#ifndef FEASST_PATCH_PATCH_ANGLE_H_
#define FEASST_PATCH_PATCH_ANGLE_H_

#include "configuration/include/model_params.h"

namespace feasst {

class PatchAngle : public ModelParam {
 public:
  PatchAngle() : ModelParam() { class_name_ = "patch_angle"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<PatchAngle>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<PatchAngle>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit PatchAngle(std::istream& istr);
  virtual ~PatchAngle() {}
};

class CosPatchAngle : public ModelParam {
 public:
  CosPatchAngle() : ModelParam() { class_name_ = "cos_patch_angle"; }
  double compute(const int type, const ModelParams& model_params) override;
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<CosPatchAngle>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<CosPatchAngle>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit CosPatchAngle(std::istream& istr);
  virtual ~CosPatchAngle() {}
};

class Director : public ModelParam {
 public:
  Director() : ModelParam() { class_name_ = "director"; }
  double compute(const int type, const ModelParams& model_params) override;
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Director>(istr); }
//  std::shared_ptr<ModelParam> create(argtype * args) const override {
//    return std::make_shared<Director>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Director(std::istream& istr);
  virtual ~Director() {}
};

}  // namespace feasst

#endif  // FEASST_PATCH_PATCH_ANGLE_H_
