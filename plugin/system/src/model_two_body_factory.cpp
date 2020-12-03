#include <memory>
#include "system/include/model_two_body_factory.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapModelTwoBodyFactory {
 public:
  MapModelTwoBodyFactory() {
    ModelTwoBodyFactory().deserialize_map()["ModelTwoBodyFactory"] =
      MakeModelTwoBodyFactory();
  }
};

static MapModelTwoBodyFactory mapper_ = MapModelTwoBodyFactory();

ModelTwoBodyFactory::ModelTwoBodyFactory(
    std::vector<std::shared_ptr<ModelTwoBody> > models)
  : ModelTwoBodyFactory() {
  for (auto model : models) {
    add(model);
  }
}

void ModelTwoBodyFactory::precompute(const ModelParams& existing) {
  for (const std::shared_ptr<Model> model : models_) {
    model->precompute(existing);
  }
}

double ModelTwoBodyFactory::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  double en = 0;
  // for (const std::shared_ptr<Model> model : models_) {
  for (int im = 0; im < num(); ++im) {
    en += models_[im]->energy(squared_distance, type1, type2, model_params);
  }
//    int index = 0;
//    while ((index < static_cast<int>(models_.size())) and
//           (en < NEAR_INFINITY)) {
//      en += models_[index]->energy(squared_distance,
//                                   type1, type2, model_params);
//      ++index;
//    }
  // DEBUG("en " << en);
  return en;
}

void ModelTwoBodyFactory::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(573, ostr);
  feasst_serialize_fstdr(models_, ostr);
}

ModelTwoBodyFactory::ModelTwoBodyFactory(std::istream& istr)
  : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(573 == version, version);
  { // HWH for unknown reasons, template deserialization fails here.
    int dim1;
    istr >> dim1;
    models_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstdr(models_[index], istr);
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

}  // namespace feasst
