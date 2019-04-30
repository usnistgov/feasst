
#ifndef FEASST_CORE_MODEL_SQUARE_WELL_H_
#define FEASST_CORE_MODEL_SQUARE_WELL_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

class ModelSquareWell : public ModelTwoBody {
 public:
  ModelSquareWell() {}

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const override {
    const double& sigma = model_params.mixed_sigma()[type1][type2];
    const double& epsilon = model_params.mixed_epsilon()[type1][type2];
    if (squared_distance <= sigma*sigma) {
      return NEAR_INFINITY;
    }
    return -epsilon;
  }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelSquareWell>(istr);
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

  ModelSquareWell(std::istream& istr) {
    feasst_deserialize_version(istr);
  }

  virtual ~ModelSquareWell() {}

 private:
  const std::string class_name_ = "ModelSquareWell";
};

inline std::shared_ptr<ModelSquareWell> MakeModelSquareWell() {
  return std::make_shared<ModelSquareWell>();
}

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_SQUARE_WELL_H_
