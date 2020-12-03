
#ifndef FEASST_MODELS_SQUARE_WELL_H_
#define FEASST_MODELS_SQUARE_WELL_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  The following model parameters are utilized:
  The sigma parameter is the hard sphere.
  The cutoff is the range of the attraction.
  The epsilon parameter is the well depth.
 */
class SquareWell : public ModelTwoBody {
 public:
  SquareWell() { class_name_ = "SquareWell"; }

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<SquareWell>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit SquareWell(std::istream& istr);
  virtual ~SquareWell() {}
};

inline std::shared_ptr<SquareWell> MakeSquareWell() {
  return std::make_shared<SquareWell>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_SQUARE_WELL_H_
