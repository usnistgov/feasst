
#ifndef FEASST_SYSTEM_MODEL_TWO_BODY_FACTORY_H_
#define FEASST_SYSTEM_MODEL_TWO_BODY_FACTORY_H_

#include "system/include/model_two_body.h"

namespace feasst {

class ModelTwoBodyFactory : public ModelTwoBody {
 public:
  ModelTwoBodyFactory() {}

  void add_model(std::shared_ptr<ModelTwoBody> model) {
    models_.push_back(model); }
  const std::vector<std::shared_ptr<Model> >& models() const {
    return models_; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    double en = 0;
    for (const std::shared_ptr<Model> model : models_) {
      en += model->energy(squared_distance, type1, type2, model_params);
    }
//    int index = 0;
//    while ((index < static_cast<int>(models_.size())) and
//           (en < NEAR_INFINITY)) {
//      en += models_[index]->energy(squared_distance,
//                                   type1, type2, model_params);
//      ++index;
//    }
    //INFO("en " << en);
    return en;
  }

  void precompute(const ModelParams& existing) override {
    for (const std::shared_ptr<Model> model : models_) {
      model->precompute(existing);
    }
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTwoBodyFactory>(istr); }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(573, ostr);
    feasst_serialize_fstdr(models_, ostr);
  }

  ModelTwoBodyFactory(std::istream& istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(573 == version, version);
    { // HWH for unknown reasons, template deserialization fails here.
      int dim1;
      istr >> dim1;
      models_.resize(dim1);
      for (int index = 0; index < dim1; ++index) {
        //feasst_deserialize_fstdr(models_[index], istr);
        {
          int existing;
          istr >> existing;
          if (existing != 0) {
            models_[index] = models_[index]->deserialize(istr);
          }
        }
      }
    }
  }

 private:
  const std::string class_name_ = "ModelTwoBodyFactory";
  std::vector<std::shared_ptr<Model> > models_;
};

inline std::shared_ptr<ModelTwoBodyFactory> MakeModelTwoBodyFactory() {
  return std::make_shared<ModelTwoBodyFactory>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_TWO_BODY_FACTORY_H_
