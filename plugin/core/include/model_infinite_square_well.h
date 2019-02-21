
#ifndef FEASST_CORE_MODEL_INFINITE_SQUARE_WELL_H_
#define FEASST_CORE_MODEL_INFINITE_SQUARE_WELL_H_

#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

/**
  U(r) = 0 when |r| < width/2, otherwise infinity.
 */
class ModelInfiniteSquareWell : public ModelTwoBody {
 public:
  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) const override {
    if (squared_distance > width_*width_/4.) {
      return NEAR_INFINITY;
    }
    return 0;
  }
  virtual ~ModelInfiniteSquareWell() {}

  void set_width(const double width) { width_ = width; }
  double width() const { return width_; }

 private:
  double width_;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_INFINITE_SQUARE_WELL_H_
