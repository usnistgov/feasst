
#ifndef FEASST_SYSTEM_MODEL_EMPTY_H_
#define FEASST_SYSTEM_MODEL_EMPTY_H_

#include "system/include/model_one_body.h"

namespace feasst {

/// This model throws an exception if used.
class ModelEmpty : public ModelOneBody {
 public:
  ModelEmpty() {}

  double energy(
      const Site& site,
      const Configuration * config,
      const ModelParams& model_params) const override {
    ERROR("Empty model should not be called");
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelEmpty>(istr);
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(189, ostr);
  }

  ModelEmpty(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(189 == version, version);
  }

  virtual ~ModelEmpty() {}

 private:
  const std::string class_name_ = "ModelEmpty";
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_EMPTY_H_
