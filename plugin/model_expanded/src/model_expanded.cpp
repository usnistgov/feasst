#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "model_expanded/include/model_expanded.h"

namespace feasst {

class MapModelExpanded {
 public:
  MapModelExpanded() {
    ModelExpanded().deserialize_map()["ModelExpanded"] =
      MakeModelExpanded();
  }
};

static MapModelExpanded mapper_ = MapModelExpanded();

ModelExpanded::ModelExpanded(argtype * args) : ModelTwoBodyFactory(args) {
  class_name_ = "ModelExpanded";
  model_index_ = integer("model_index", args, 0);
}
ModelExpanded::ModelExpanded(argtype args) : ModelExpanded(&args) {
  feasst_check_all_used(args);
}

double ModelExpanded::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  DEBUG("model_index " << model_index_);
  DEBUG("num models " << num());
  ASSERT(model_index_ >= 0, "model_index:" << model_index_ << " < 0");
  ASSERT(model_index_ < num(),
    "model_index:" << model_index_ << " >= number of models:" << num());
  return models_[model_index_]->energy(squared_distance, type1, type2, model_params);
}

void ModelExpanded::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_two_body_factory_(ostr);
  feasst_serialize_version(7168, ostr);
  feasst_serialize(model_index_, ostr);
}

ModelExpanded::ModelExpanded(std::istream& istr)
  : ModelTwoBodyFactory(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(7168 == version, version);
  feasst_deserialize(&model_index_, istr);
}

void ModelExpanded::set_model_index(const int index) {
  if (index >= 0 && index < static_cast<int>(num())) {
    model_index_ = index;
  }
}

}  // namespace feasst
