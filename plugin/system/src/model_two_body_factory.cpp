#include <memory>
#include "utils/include/arguments.h"
#include "utils/include/file.h"
#include "utils/include/serialize.h"
#include "system/include/model_two_body_factory.h"

namespace feasst {

class MapModelTwoBodyFactory {
 public:
  MapModelTwoBodyFactory() {
    ModelTwoBodyFactory().deserialize_map()["ModelTwoBodyFactory"] =
      MakeModelTwoBodyFactory();
  }
};

static MapModelTwoBodyFactory mapper_ = MapModelTwoBodyFactory();

ModelTwoBodyFactory::ModelTwoBodyFactory(argtype * args) {
  class_name_ = "ModelTwoBodyFactory";
  std::string model_file = str("model_file", args, "");
  if (!model_file.empty()) {
    ASSERT(!used("model0", *args),
      "ModelTwoBodyFactory should not be given both model_file and model0");
    ASSERT(file_exists(model_file),
      "model_file: " << model_file << "not found");
    std::ifstream file(model_file);
    ASSERT(file.good(), "error");
    const bool is_found = find("ModelTwoBodyFactory", file);
    ASSERT(is_found, "ModelTwoBodyFactory not found in " << model_file);
    std::string line;
    std::getline(file, line);
    ASSERT(line.empty(), "error");
    while (std::getline(file, line)) {
      DEBUG(line);
      std::pair<std::string, argtype> margs = parse_line(line, NULL, NULL);
      auto model = ModelTwoBody().factory(margs.first, &margs.second);
      FEASST_CHECK_ALL_USED(margs.second);
      models_.push_back(model);
    }
  } else {
    int model_index = 0;
    std::stringstream key;
    key << "model" << model_index;
    while (used(key.str(), *args)) {
      const std::string model_name = str(key.str(), args);
      auto model = ModelTwoBody().factory(model_name, args);
      models_.push_back(model);
      ++model_index;
      key.str("");
      key << "model" << model_index;
    }
  }
  DEBUG("num " << num());
}
ModelTwoBodyFactory::ModelTwoBodyFactory(argtype args) : ModelTwoBodyFactory(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void ModelTwoBodyFactory::add(
    std::vector<std::shared_ptr<ModelTwoBody> > models) {
  for (auto model : models) {
    add(model);
  }
}

void ModelTwoBodyFactory::precompute(const ModelParams& existing) {
  for (const std::shared_ptr<Model>& model : models_) {
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
  serialize_model_two_body_factory_(ostr);
}

void ModelTwoBodyFactory::serialize_model_two_body_factory_(std::ostream& ostr) const {
  serialize_model_(ostr);
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

std::shared_ptr<ModelTwoBodyFactory> MakeModelTwoBodyFactory(
    std::shared_ptr<ModelTwoBody> model1,
    std::shared_ptr<ModelTwoBody> model2) {
  auto fac = std::make_shared<ModelTwoBodyFactory>();
  fac->add(model1);
  fac->add(model2);
  return fac;
}

}  // namespace feasst
