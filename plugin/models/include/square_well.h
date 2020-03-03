
#ifndef FEASST_MODELS_SQUARE_WELL_H_
#define FEASST_MODELS_SQUARE_WELL_H_

#include "system/include/model_two_body.h"
#include "math/include/constants.h"

namespace feasst {

/**
  The following model parameters are utilized:
  The sigma parameter is the hard sphere.
  The cutoff is the range of the attraction.
  The epsilon parameter is the well depth.
 */
class SquareWell : public ModelTwoBody {
 public:
  SquareWell() {}

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
    return std::make_shared<SquareWell>(istr);
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
  }

  SquareWell(std::istream& istr) {
    feasst_deserialize_version(istr);
  }

  virtual ~SquareWell() {}

 private:
  const std::string class_name_ = "SquareWell";
};

inline std::shared_ptr<SquareWell> MakeSquareWell() {
  return std::make_shared<SquareWell>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_SQUARE_WELL_H_
